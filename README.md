# linearized-NSE-data-assimilation

# How to reproduce 
Parameters to be changed by the user are located in the fourth block titled as "setting parameters":
* `order` describes the polynomial degree of the finite elements 
* `domain_case` which describes the geometrical setup (convex, non-convex ..)  
* `theta` describes the perturbation order (as defined in the paper) 
* `give_pressure` is a flag which describes whether global pressure data is added tot the Lagrangian (see Section 5.2)

## Figure 2 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = None` 
* `domain_case = subdomain_cases[2]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 2a can be found in the folder `plots` and is called `stokes-convex-k123.tex`.

## Figure 2 b) 

