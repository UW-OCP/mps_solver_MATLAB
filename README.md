# mps_solver_MATLAB

Multiple shooting method for solving optimal control problems.

The repo contains two directories: 1)ms_BVP; 2) ms_OCP.

```html
<img src="http://chart.googleapis.com/chart?cht=tx&chl=\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" style="border:none;">
- <img src="https://latex.codecogs.com/gif.latex?O_t=\text { Onset event at time bin } t " />
<img src="https://latex.codecogs.com/gif.latex?O_t=\text { Onset event at time bin } t " />
```

## mps_BVP

The directory contains the solver that can solves the boundary value problem (BVP) with differential algebraic equations (DAEs) defined in the specific format and gives the solution of the BVP.

The input of the solver should contain the various function handles to define the BVP and the initial estimate to the problem. The solver gives the solution to the BVP.

### Input

* ODE_h : function handle to define the ordinary differential equations (ODEs). The API of the function should be $$\dot{y} = ODE\_h(y, z, p, \alpha)$$. The dynamic systems considered here are non-autonomous system so there is no explicitly time $t$ in the ODEs. The variables of the function are defined as below,
  * Input:
     * $y$ : ODE variables.
     * $z$ : DAE variables.
     * $p$ : parameter variables.
     * $\alpha$ : continuation parameter.
  * Output:
    * $\dot{y}$ : derivatives of the ODE variables.
* DAE_g : function handle to define the differential algebraic equations (DAEs). The API of the function should be $0 = DAE\_g(y, z, p, \alpha)$. The variables of the function are defined as below, 
  * Input:
    * $y$ : ODE variables.
    * $z$ : DAE variables.
    * $p$ : parameter variables.
    * $\alpha$ : continuation parameter.
* BC_r : function handle to define the boundary constraints of the BVP problem. The API of the function should be $0 = BC{\_}r(y_{0}, y_{M}, p)$. The variables of the function are defined as below,
  * Input:
    * $y_{0}$ : values of the ODE variables at the initial time $t_{0}$.
    * $y_{M}$ : values of the ODE variables at the final time $t_{M}$.
    * $p​$ : parameter variables.
* D_hg : function handle to define the derivative of the ODEs and DAEs w.r.t. ODE variables, DAE variables, and parameter variables. The API of the function should be $[h_{y}, h_{z}, h_{p}, g_{y}, g_{z}, g_{p}] = D{\_}hg(s, p, \alpha)​$. The variables of the function are defined as below,
  * Input:
    * $s$ : vector form of ODE variables and DAE variables.
    * $p$ : parameter variables.
    * $\alpha$ : continuation parameter.
  * Output:
    * $h_y​$ : derivative of the ODEs $h​$ w.r.t. the ODE variables $y​$.
    * $h_z$ : derivative of the ODEs $h$ w.r.t. the ODE variables $z$.
    * $h_p$ : derivative of the ODEs $h$ w.r.t. the ODE variables $p$.
    * $g_y$ : derivative of the DAEs $g$ w.r.t. the ODE variables $y$.
    * $g_z$ : derivative of the DAEs $g$ w.r.t. the ODE variables $z$.
    * $g_p$ : derivative of the DAEs $g$ w.r.t. the ODE variables $p$.
* D_r : function handle to define the derivative of the boundary constraints w.r.t. the initial and final ODE variable values. The API of the function should be $[r_{y_{0}}, r_{y_{M}}, r_{p}] = D{\_}r(y_{0}, y_{M}, p)​$. The variables of the function are defined as below,
  * Input:
    * $y_{0}$ : values of the ODE variables at the initial time $t_{0}$.
    * $y_{M}$ : values of the ODE variables at the final time $t_{M}$.
    * $p$ : parameter variables.
  * Output:
    * $r_{y_{0}}$ : derivative of the boundary constraints w.r.t. the initial ODE variable values.
    * $r_{y_{M}}$ : derivative of the boundary constraints w.r.t. the final ODE variable values.
    * $r_{p}$ : derivative of the boundary constraints w.r.t. the parameter variables.
* y0 : matrix form of the intial estimate of the ODE variables $y$ at each time node, where each row corresponds to the values at the same node, the dimension of the matrix should be $N \times n_{y}$, where $N$ is the number of time nodes and $n_{y}$ is the number of ODE variables.
*  z0​ :  matrix form of the intial estimate of the DAE variables $z$ at each time node, where each row corresponds to the values at the same node, the dimension of the matrix should be $N \times n_{z}$, where $N$ is the number of time nodes and $n_{z}$ is the number of DAE variables.
* p​ : vector form of the initial estimate of the parameter variables, the dimension of the vector should be $n_{p}​$, where $n_{p}​$ is the number of parameter variables.
* ti : initial time of the problem, which should be a nonnegative number.
* tM : final time of the problem, which should be a positive number.
* output : file name of the output files which should be a string.
* alpham : the termination tolerance of the continuation parameter.
* tolerance : numerical tolerance of the numerical iterations.
* nodes_max : number of maximum nodes allowed during the iterations.
* mesh_max : number of maximum mesh refinement allowed during the iterations.
* nx : number of state variables of the problem.
* nu : number of control variables of the problem.
* nw : number of parameters of the problem.
* nd : number of control variable inequality constraints (CVICs) of the problem.
* ns : number of state variable inequality constraints (SVICs) of the problem.

### Output

* y : final solution of the ODE varaibles of the problem which is saved in output + '_y.txt' file.
* z : final solution of the DAE varaibles of the problem which is saved in output + '_z.txt' file.
* p : final solution of the parameter varaibles of the problem which is saved in output + '_p.txt' file.
* $\alpha$ : final value of the continuation parameter which in saved in output + '_alpha.txt' file.

## mps_OCP

The directory contains the solver that can directly solves the optimal control problem (OCP) defined in the specific format and gives the optimal solution.

The input is just the plain text file with the specific defined name fileds and the ouputs are the optimal solutions of various variables of the OCP.