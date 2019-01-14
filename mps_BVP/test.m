clear variables;

N = 101;
y0 = ones(N, 4);
z0 = ones(N, 1);
para = ones(4, 1);
t0 = 0;
tf = 1;
output_File = 'ex2data';
alpham = 1;
tol = 1e-6;
nodes_max = 2000;
mesh_max = 25;
nx = 2;
nu = 1;
nw = 0;
nd = 0;
ns = 0;

mps_BVP(@(y, z, p, alpha) ODE_hh(y, z, p, alpha), @(y, z, p, alpha) DAE_gg(y, z, p, alpha),...
    @(y0, yM, p) boundary_constraint(y0, yM, p), @(s, p, alpha) difference_DAE(s, p, alpha),...
    @(y0, yM, p) difference_BC(y0, yM, p), @(t, s, p, alpha) jacobian_DAE(t, s, p, alpha),...
    y0, z0, para, t0, tf, output_File, alpham, tol, nodes_max, mesh_max, nx, nu, nw, nd, ns);