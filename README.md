# linearized-NSE-data-assimilation

# README 
This repository contains the software, data and instructions to reproduce the numerical experiments in the paper 
> Data assimilation finite element method for the linearized Navier-Stokes equations with higher order polynomial approximation
> 
> * authors: Erik Burman, Deepika Garg, Janosch Preuss
> * University College London 

# How to run / install
* Please download the docker image from zenodo.
* Assuming that `linearized-nse-repro.tar` is the filename of the downloaded image, please load the image with `docker load < linearized-nse-repro.tar`. Usually this command has to be executed with root privileges, which means e.g. on linux systems that `sudo` has to be added in front of this command.
* Run the image with `docker run --init -ti -p 8888:8888 linearized-nse-repro:v1`. Open the url shown in the terminal in your browser.
* Proceed further as described in [How to reproduce](#repro)

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
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 2a can be 
found in the folder `plots` and is called `stokes-convex-k123.tex`.

## Figure 2 b) Run the same file as Figure 2 a) the only changes are:
* `domain_case = subdomain_cases[3]`

* The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 2b can be 
found in the folder `plots` and is called `stokes-nonconvex-k123.tex`.


## Figure 3 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 0` 
* `domain_case = subdomain_cases[2]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 3a can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation0.tex`.


## Figure 3 b) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 1` 
* `domain_case = subdomain_cases[2]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 3b can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation1.tex`.

## Figure 3 c) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 2` 
* `domain_case = subdomain_cases[2]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 3c can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation2.tex`.


## Figure 4 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 0` 
* `domain_case = subdomain_cases[3]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 4a can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation0.tex`.


## Figure 4 b) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 1` 
* `domain_case = subdomain_cases[3]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 4b can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation1.tex`.

## Figure 4 c) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 2` 
* `domain_case = subdomain_cases[3]` 
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 4a can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation2.tex`.


## Figure 5 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = None` 
* `domain_case = subdomain_cases[2]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 5a can be 
found in the folder `plots` and is called `stokes-convex-k123.tex`.

## Figure 5 b) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = None` 
* `domain_case = subdomain_cases[3]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 5b can be 
found in the folder `plots` and is called `stokes-nonconvex-k123.tex`.


## Figure 6 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 0` 
* `domain_case = subdomain_cases[2]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 6a can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation0.tex`.

## Figure 6 b) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 1` 
* `domain_case = subdomain_cases[2]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 6b can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation1.tex`.

## Figure 6 c) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 0` 
* `domain_case = subdomain_cases[2]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-convex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 6c can be 
found in the folder `plots` and is called `stokes-convex-k123-data-perturbation2.tex`.


## Figure 7 a) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 0` 
* `domain_case = subdomain_cases[3]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 7a can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation0.tex`.



## Figure 7 b) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 1` 
* `domain_case = subdomain_cases[3]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 7b can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation1.tex`.

## Figure 7 c) 
Run the file `Ill_posed_Stokes.ipynb` using the following parameters:  
* `theta = 2` 
* `domain_case = subdomain_cases[3]` 
* `dual_finite_element_u  = ufl.VectorElement("CG", mesh.ufl_cell(), 1)`
* `primal_finite_element_p = ufl.FiniteElement("CG", mesh.ufl_cell(), order-1)`
* `dual_finite_element_p  = ufl.FiniteElement("CG", mesh.ufl_cell(), 1)`
For each polynomial degree `order` in `[2,3]` a separate run is required. 

The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
stokes-ill-posed-nonconvex-Helmholtz-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 7c can be 
found in the folder `plots` and is called `stokes-nonconvex-k123-data-perturbation2.tex`.


## Figure 9 a) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = None` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  1`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 9a can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123.tex`.


## Figure 9 b) Run the same file as Figure 9 a) only change is:

* `give_pressure = True` 

## Figure 9 c) Run the same file as Figure 9 a) only change is:

* ` diffusion_coff=  0.01`

## Figure 9 d) Run the same file as Figure 9 c) only change is:

* `give_pressure = True` 

## Figure 9 e) Run the same file as Figure 9 a) only change is:

* ` diffusion_coff=  0.0001`

## Figure 9 f) Run the same file as Figure 9 e) only change is:

* `give_pressure = True` 

## Figure 9 g) Run the same file as Figure 9 a) only change is:

* ` diffusion_coff=  0`

## Figure 9 h) Run the same file as Figure 9 g) only change is:

* `give_pressure = True` 


## Figure 10 a) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 0` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  1`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 10a can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation0.tex`.


## Figure 10 b) Run the same file as Figure 10 a) only change is:

* `give_pressure = True` 

## Figure 10 c) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 0` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  0.01`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 10c can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation0.tex`.


## Figure 10 d) Run the same file as Figure 10 c) only change is:

* `give_pressure = True` 

## Figure 10 e) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 0` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  0.0001`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 10e can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation0.tex`.


## Figure 10 f) Run the same file as Figure 10 e) only change is:

* `give_pressure = True` 


## Figure 10 g) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 0` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  0`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 10g can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation0.tex`.


## Figure 10 h) Run the same file as Figure 10 g) only change is:

* `give_pressure = True` 


## Figure 11 a) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 1` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  1`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 11a can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation1.tex`.


## Figure 11 b) Run the same file as Figure 11 a) only change is:

* `give_pressure = True` 


## Figure 11 c) 
Run the file `Ill_posed_new_Navier_Stokes` using the following parameters:  
* `theta = 2` 
* `give_pressure = False` 
* `domain_case = subdomain_cases[5]` 
* ` diffusion_coff=  1`
For each polynomial degree `order` in `[1,2,3]` a separate run is required. 
The convergence can be seen in the plot generated as the end of the jupyter notebook.
The data for generating the plot included in the paper will be saved to the folder `data` in `.dat` files named `
Navier-Stokes-ill-posed-flow_downstream-order__j__-theta__i__.dat` where __j__ and __i__ represent the corresponding 
values of `theta` and `order respectively`. The `.tex` file for generating the Fig. 11c can be 
found in the folder `plots` and is called `Navier-stokes-downstream-k123-data-perturbation2.tex`.


## Figure 11 d) Run the same file as Figure 11 c) only change is:

* `give_pressure = True` 
