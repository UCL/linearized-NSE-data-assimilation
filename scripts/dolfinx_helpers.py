from typing import Callable, Optional, Union
import numpy as np
import gmsh
import pygmsh
import ufl
from dolfinx.fem import (
    Function,
    FunctionSpace,
    Expression,
    dirichletbc,
    DirichletBCMetaClass,
    Constant,
    locate_dofs_topological,
)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import (
    extract_gmsh_geometry,
    extract_gmsh_topology_and_markers,
    ufl_mesh_from_gmsh,
)
from dolfinx.mesh import create_mesh, compute_boundary_facets
from dolfinx.cpp.io import perm_gmsh
from dolfinx.cpp.mesh import to_type, Mesh
from mpi4py import MPI
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import Axes3D


def create_dolfinx_mesh_using_pygmsh(
    define_geometry: Callable[[pygmsh.common.CommonGeometry], None],
    dim: int = 2,
    geometry_class: pygmsh.common.CommonGeometry = pygmsh.geo.Geometry,
    **kwargs
) -> Mesh:
    """Create a dolfinx mesh object using pygmsh.

    Args:
        define_geometry: Function which given a pygmsh geometry instance,
            defines the mesh geometry using the object's methods.
        dim: Spatial dimension.
        geometry_class: pygmsh geometry class to use.
        kwargs: Any keyword arguments to pass to `define_geometry`.

    Returns:
        dolfinx finite element mesh with specified geometry.
    """
    with geometry_class() as geom:
        define_geometry(geom, **kwargs)
        # Generate mesh using pygmsh - this returns a meshio mesh object
        # but also updates the gmsh.model global object which we will
        # use directly here to allow using the helper functions in
        # dolfinx.io to directly load this mesh into dolfinx
        geom.generate_mesh(dim=dim)
        # Get mesh geometry
        geometry_data = extract_gmsh_geometry(gmsh.model)
        # Get mesh topology for each element
        topology_data = extract_gmsh_topology_and_markers(gmsh.model)
        # Extract the cell type and number of nodes per cell
        cell_id = next(iter(topology_data.keys()))
        properties = gmsh.model.mesh.getElementProperties(cell_id)
    name, dim, order, num_nodes, local_coords, _ = properties
    cells = topology_data[cell_id]["topology"]
    # Permute topology data from MSH-ordering to dolfinx-ordering
    ufl_domain = ufl_mesh_from_gmsh(cell_id, dim)
    gmsh_cell_perm = perm_gmsh(to_type(str(ufl_domain.ufl_cell())), num_nodes)
    cells = cells[:, gmsh_cell_perm]
    # Create dolfinx mesh
    return create_mesh(MPI.COMM_WORLD, cells, geometry_data[:, :dim], ufl_domain)


def get_matplotlib_triangulation_from_mesh(mesh: Mesh) -> Triangulation:
    """Get matplotlib triangulation corresponding to dolfinx mesh.

    Args:
        mesh: Finite element mesh to get triangulation for.

    Returns:
        Object representing triangulation of mesh to use in Matplotlib plot functions.
    """
    assert mesh.topology.dim == 2, "Only two-dimensional spatial domains are supported"
    # The triangulation of the mesh corresponds to the connectivity between elements of
    # dimension 2 (triangles) and elements of dimension 0 (points)
    mesh.topology.create_connectivity(2, 0)
    triangles = mesh.topology.connectivity(2, 0).array.reshape((-1, 3))
    return Triangulation(
        mesh.geometry.x[:, 0],
        mesh.geometry.x[:, 1],
        triangles,
    )


def project_expression_on_function_space(
    expression: ufl.core.expr.Expr, function_space: ufl.FunctionSpace
) -> Function:
    """Project expression onto finite element function space.

    Args:
        expression: UFL object defining expression to project.
        function_space: Finite element function space.

    Returns:
        Function representing projection of expression.
    """
    trial_function = ufl.TrialFunction(function_space)
    test_function = ufl.TestFunction(function_space)
    return LinearProblem(
        ufl.inner(trial_function, test_function) * ufl.dx,
        ufl.inner(expression, test_function) * ufl.dx,
    ).solve()


def plot_functions_as_heatmaps(
    mesh: Mesh,
    functions_dict: dict[str, Function],
    ax_size: float = 5.0,
    show_colorbar: bool = True,
):
    """Plot one or more finite element functions on 2D domains as heatmaps.

    Args:
        mesh: Finite element mesh corresponding to domain of function(s).
        functions_dict: Dictionary from string labels to `dolfinx` functions.
        ax_size: Size of axis to plot each function on in inches.
        show_colorbar: Whether to show colorbar key for function values.
        show_triangulation: Whether to show triangulation.

    Return:
        Tuple with first entry Matplotlib figure and second (array of) axes.
    """
    assert mesh.topology.dim == 2, "Only two-dimensional spatial domains are supported"
    num_func = len(functions_dict)
    multiplier = 1.25 if show_colorbar else 1.0
    fig, axes = plt.subplots(
        1, num_func, figsize=(multiplier * num_func * ax_size, ax_size)
    )
    triangulation = get_matplotlib_triangulation_from_mesh(mesh)
    for ax, (label, func) in zip(np.atleast_1d(axes), functions_dict.items()):
        artist = ax.tripcolor(
            triangulation,
            func.x.array.real,
        )
        ax.set_title(label)
        ax.set_aspect(1)
        ax.set(xlabel="Spatial coordinate 0", ylabel="Spatial coordinate 1")
        if show_colorbar:
            fig.colorbar(artist, ax=ax, fraction=0.125, pad=0.05)
    fig.tight_layout()
    return fig, axes


def plot_functions_as_surfaces(
    mesh: Mesh, functions_dict: dict[str, Function], ax_size: float = 5.0
):
    """Plot one or more finite element functions on 2D domains as surfaces.

    Args:
        mesh: Finite element mesh corresponding to domain of function(s).
        functions_dict: Dictionary from string labels to `dolfinx` functions.
        ax_size: Size of axis to plot each function on in inches.

    Return:
        Tuple with first entry Matplotlib figure and second (array of) axes.
    """
    assert mesh.topology.dim == 2, "Only two-dimensional spatial domains are supported"
    num_func = len(functions_dict)
    fig, axes = plt.subplots(
        1,
        num_func,
        figsize=(num_func * ax_size, ax_size),
        subplot_kw={"projection": "3d"},
    )
    triangulation = get_matplotlib_triangulation_from_mesh(mesh)
    for ax, (label, func) in zip(np.atleast_1d(axes), functions_dict.items()):
        ax.plot_trisurf(
            triangulation,
            func.compute_point_values()[:, 0],
        )
        ax.set(
            xlabel="Spatial coordinate 0", ylabel="Spatial coordinate 1", zlabel=label
        )
    fig.tight_layout()
    return fig, axes


def define_dirchlet_boundary_conditions(
    boundary_values: Union[Function, Constant], 
    function_space: Optional[FunctionSpace] = None
) -> DirichletBCMetaClass:
    """Define dolfinx object representing Dirichlet boundary conditions.
    
    Args:
        boundary_values: dolfinx function specifying values across
            domain from which fixed values at boundary will be taken.
        function_space: Optional argument specifying finite element
            function space from which boundary degrees of freedom
            should be computed. If `None` (default) then
            `boundary_values.function_space` is used (which will
            only work if `boundary_values` is a `Function` instance.

    Returns:
        Dirichlet boundary condition object.
    """ 
    if function_space is None:
        function_space = boundary_values.function_space
    mesh = function_space.mesh
    # dolfinx no longer allows passing a constant directly
    # if not isinstance(boundary_values, Function):
    #     boundary_function = Function(function_space)
    #     boundary_function.x.array[:] = boundary_values
    #     boundary_values = boundary_function
    # Create facet to cell connectivity required to determine boundary facets
    facet_dim = mesh.topology.dim - 1
    mesh.topology.create_connectivity(facet_dim, mesh.topology.dim)
    boundary_facets = np.flatnonzero(compute_boundary_facets(mesh.topology))
    boundary_dofs = locate_dofs_topological(function_space, facet_dim, boundary_facets)
    return dirichletbc(boundary_values, boundary_dofs, function_space)


def define_poisson_equation_variational_forms(
    trial_function: ufl.Argument,
    test_function: ufl.Argument,
    source_term: ufl.core.expr.Expr,
) -> tuple[ufl.Form, ufl.Form]:
    """Define the bilinear and linear forms for the variational form of Poisson equation

    \[ \int_{\Omega} \nabla u \dot \nabla v dx = \int_{\Omega} f v dx

    Args:
        trial_function: UFL object corresponding to trial function $u$.
        test_function: UFL object corresponding to test function $v$.
        source_term: UFL object corresponding to source term $f$.

    Returns:
        Tuple with first entry corresponding to bilinear form and second entry linear form.
    """
    bilinear_form = ufl.dot(ufl.grad(trial_function), ufl.grad(test_function)) * ufl.dx
    linear_form = source_term * test_function * ufl.dx
    return bilinear_form, linear_form


def define_exponential_squared_expression(
    mesh: Mesh,
    length_scale: float = 0.1,
    amplitude: float = 4.0,
    centre: tuple[float, float] = (0.0, 0.0),
) -> Expression:
    """Define expression on mesh representing exponential squared bump function

    \[ f(x) = a \exp(-((x_0 - c_0)^2 + (x_1 - c_1)^2) / s^2) \]

    where $a$ is an amplitude parameter, $c$ a 2-tuple specifying the centre
    point and $s$ a length scale parameter.

    Args:
        mesh: Finite element mesh to define expression on.
        length_scale: Parameter defining extent of bump in spatial dimension.
        amplitude: Parameter defining extent of bump in 'function' dimension.
        centre: Location of centre of bump in spatial domain.

    Return:
        Expression representing bump function on domain.
    """
    x = ufl.SpatialCoordinate(mesh)
    assert len(x) == len(centre), "Dimension of centre does not match mesh"
    amplitude = Constant(mesh, amplitude)
    centre = Constant(mesh, centre)
    return amplitude * ufl.exp(
        -sum((x[i] - centre[i]) ** 2 for i in range(mesh.topology.dim))
        / length_scale ** 2
    )


def project_expression_on_function_space(
    expression: ufl.core.expr.Expr, function_space: ufl.FunctionSpace
) -> Function:
    """Project expression onto finite element function space.

    Args:
        expression: UFL object defining expression to project.
        function_space: Finite element function space.

    Returns:
        Function representing projection of expression.
    """
    trial_function = ufl.TrialFunction(function_space)
    test_function = ufl.TestFunction(function_space)
    return LinearProblem(
        ufl.inner(trial_function, test_function) * ufl.dx,
        ufl.inner(expression, test_function) * ufl.dx,
    ).solve()


def solve_linear_problem(
    bilinear_form: ufl.Form,
    linear_form: ufl.Form,
    boundary_conditions: list[DirichletBCMetaClass],
    **kwargs
) -> Function:
    """Solve a linear variation problem.

    Args:
        bilinear_form: UFL object representing bilinear form of problem.
        linear_form: UFL object representing linear form of problem.
        boundary_conditions: Any Dirichlet boundary conditions in problem.
        kwargs: Any additional keyword arguments to pass to `LinearProblem` initialiser.

    Returns:
        Finite element function corresponding to solution to problem.
    """
    return LinearProblem(
        bilinear_form, linear_form, bcs=boundary_conditions, **kwargs
    ).solve()
