# mps_solver_MATLAB
Multiple shooting method for solving optimal control problems.

The repo contains two directories: 1)ms_BVP; 2) ms_OCP.

## mps_BVP

The directory contains the solver that can solves the boundary value problem (BVP) defined in the specific format and gives the solution of the BVP.

The input of the solver should contain the various function handles to define the BVP and the initial estimate to the problem. The solver gives the solution to the BVP.

## mps_OCP

The directory contains the solver that can directly solves the optimal control problem (OCP) defined in the specific format and gives the optimal solution.

The input is just the plain text file with the specific defined name fileds and the ouputs are the optimal solutions of various variables of the OCP.