# mps_solver_MATLAB

Multiple shooting method for solving optimal control problems.

The repo contains two directories: 1)ms_BVP; 2) ms_OCP.

## mps_BVP

The directory contains the solver that can solves the boundary value problem (BVP) with differential algebraic equations (DAEs) defined in the specific format and gives the solution of the BVP.

The input of the solver should contain the various function handles to define the BVP and the initial estimate to the problem. The solver gives the solution to the BVP.

### Input

* ODE_h : function handle to define the ordinary differential equations (ODEs). The API of the function should be ![equation](https://latex.codecogs.com/gif.latex?\dot{y}&space;=&space;ODE{\_}h(y,&space;z,&space;p,&space;\alpha)). The dynamic systems considered here are non-autonomous system so there is no explicitly time ![equation](https://latex.codecogs.com/gif.latex?t) in the ODEs. The variables of the function are defined as below,
  * Input:
     * ![equation](https://latex.codecogs.com/gif.latex?y) : ODE variables.
     * ![equation](https://latex.codecogs.com/gif.latex?z) : DAE variables.
     * ![equation](https://latex.codecogs.com/gif.latex?p) : parameter variables.
     * ![equation](https://latex.codecogs.com/gif.latex?\alpha) : continuation parameter.
  * Output:
    * ![equation](https://latex.codecogs.com/gif.latex?\dot{y}) : derivatives of the ODE variables.
* DAE_g : function handle to define the differential algebraic equations (DAEs). The API of the function should be ![equation](https://latex.codecogs.com/gif.latex?0&space;=&space;DAE{\_}g(y,&space;z,&space;p,&space;\alpha)). The variables of the function are defined as below, 
  * Input:
    * ![equation](https://latex.codecogs.com/gif.latex?y) : ODE variables.
    * ![equation](https://latex.codecogs.com/gif.latex?z) : DAE variables.
    * ![equation](https://latex.codecogs.com/gif.latex?p) : parameter variables.
    * ![equation](https://latex.codecogs.com/gif.latex?\alpha) : continuation parameter.
* BC_r : function handle to define the boundary constraints of the BVP problem. The API of the function should be ![equation](https://latex.codecogs.com/gif.latex?0&space;=&space;BC{\_}r(y_{0},&space;y_{M},&space;p)). The variables of the function are defined as below,
  * Input:
    * ![equation](https://latex.codecogs.com/gif.latex?y_{0}) : values of the ODE variables at the initial time ![equation](https://latex.codecogs.com/gif.latex?t_{0}).
    * ![equation](https://latex.codecogs.com/gif.latex?y_{M}): values of the ODE variables at the final time ![equation](https://latex.codecogs.com/gif.latex?t_{M}).
    * ![equation](https://latex.codecogs.com/gif.latex?p): parameter variables.
* D_hg : function handle to define the derivative of the ODEs and DAEs w.r.t. ODE variables, DAE variables, and parameter variables. The API of the function should be ![equation](https://latex.codecogs.com/gif.latex?[h_{y},&space;h_{z},&space;h_{p},&space;g_{y},&space;g_{z},&space;g_{p}]&space;=&space;D{\_}hg(s,&space;p,&space;\alpha)). The variables of the function are defined as below,
  * Input:
    * ![equation](https://latex.codecogs.com/gif.latex?s) : vector form of ODE variables and DAE variables.
    * ![equation](https://latex.codecogs.com/gif.latex?p) : parameter variables.
    * ![equation](https://latex.codecogs.com/gif.latex?\alpha) : continuation parameter.
  * Output:
    * ![equation](https://latex.codecogs.com/gif.latex?h_y) : derivative of the ODEs ![equation](https://latex.codecogs.com/gif.latex?h) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?y).
    * ![equation](https://latex.codecogs.com/gif.latex?h_z) : derivative of the ODEs ![equation](https://latex.codecogs.com/gif.latex?h) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?z).
    * ![equation](https://latex.codecogs.com/gif.latex?h_p) : derivative of the ODEs ![equation](https://latex.codecogs.com/gif.latex?h) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?p).
    * ![equation](https://latex.codecogs.com/gif.latex?g_y) : derivative of the DAEs ![equation](https://latex.codecogs.com/gif.latex?g) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?y).
    * ![equation](https://latex.codecogs.com/gif.latex?g_z) : derivative of the DAEs ![equation](https://latex.codecogs.com/gif.latex?g) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?z).
    * ![equation](https://latex.codecogs.com/gif.latex?g_p) : derivative of the DAEs ![equation](https://latex.codecogs.com/gif.latex?g) w.r.t. the ODE variables ![equation](https://latex.codecogs.com/gif.latex?p).
* D_r : function handle to define the derivative of the boundary constraints w.r.t. the initial and final ODE variable values. The API of the function should be ![equation](https://latex.codecogs.com/gif.latex?[r_{y_{0}},&space;r_{y_{M}},&space;r_{p}]&space;=&space;D{\_}r(y_{0},&space;y_{M},&space;p)). The variables of the function are defined as below,
  * Input:
    * ![equation](https://latex.codecogs.com/gif.latex?y_{0}) : values of the ODE variables at the initial time ![equation](https://latex.codecogs.com/gif.latex?t_{0}).
    * ![equation](https://latex.codecogs.com/gif.latex?y_{M}) : values of the ODE variables at the final time ![equation](https://latex.codecogs.com/gif.latex?t_{M}).
    * ![equation](https://latex.codecogs.com/gif.latex?p) : parameter variables.
  * Output:
    * ![equation](https://latex.codecogs.com/gif.latex?r_{y_{0}}) : derivative of the boundary constraints w.r.t. the initial ODE variable values.
    * ![equation](https://latex.codecogs.com/gif.latex?r_{y_{M}}) : derivative of the boundary constraints w.r.t. the final ODE variable values.
    * ![equation](https://latex.codecogs.com/gif.latex?r_{p}) : derivative of the boundary constraints w.r.t. the parameter variables.
* y0 : matrix form of the intial estimate of the ODE variables ![equation](https://latex.codecogs.com/gif.latex?y) at each time node, where each row corresponds to the values at the same node, the dimension of the matrix should be ![equation](https://latex.codecogs.com/gif.latex?N&space;\times&space;n_{y}), where ![equation](https://latex.codecogs.com/gif.latex?N) is the number of time nodes and ![equation](https://latex.codecogs.com/gif.latex?n_{y}) is the number of ODE variables.
*  z0 :  matrix form of the intial estimate of the DAE variables ![equation](https://latex.codecogs.com/gif.latex?z) at each time node, where each row corresponds to the values at the same node, the dimension of the matrix should be ![equation](https://latex.codecogs.com/gif.latex?N&space;\times&space;n_{z}), where ![equation](https://latex.codecogs.com/gif.latex?N) is the number of time nodes and ![equation](https://latex.codecogs.com/gif.latex?n_{z}) is the number of DAE variables.
* p : vector form of the initial estimate of the parameter variables, the dimension of the vector should be ![equation](https://latex.codecogs.com/gif.latex?n_{p}), where ![equation](https://latex.codecogs.com/gif.latex?n_{p}) is the number of parameter variables.
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
* ![equation](https://latex.codecogs.com/gif.latex?\alpha) : final value of the continuation parameter which in saved in output + '_alpha.txt' file.

## mps_OCP

The directory contains the solver that can directly solves the optimal control problem (OCP) defined in the specific format and gives the optimal solution.

The input is just the plain text file with the specific defined name fileds and the ouputs are the optimal solutions of various variables of the OCP.