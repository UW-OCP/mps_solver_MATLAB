function [y0, z0, para, t0, tf, t_span, output_File, alpham, tol, nodes_max, mesh_max, nx, nu, nw, nd, ns] = generate_bvpdae(filename)
% avoid the use of the ancient Greek alphabet like "beta", "theta"
% as the name of the variables as they are easy to be the name of matlan mupad function

%% Pre-allocate size of each variables
nx = 0;
nw = 0;
nu = 0;
nd = 0;
ns = 0;
nT = 0;
nF = 0;
nau = 0;
ntp = 0;

t0 = 0;
tf = 1;

flag_se = 0; % flag to detect whether state estimate exists in the file
flag_ce = 0; % flag to detect whether control estimate exists in the file
flag_pe = 0; % flag to detect whether parameter estimate exists in the file
flag_constants = 0; % flag to detect whether constants exists in the file
flag_inputFile = 0; % flag to detect whether inputFile exists in the file

Nodes = 101;
nodes_max = 2000;
mesh_max = 25;
alpham = 1;
output_File = 'output';

% Set up the variable to hold each term
sv = 0;
cv = 0;
ic = 0;
tc = 0;
cf = 0;
de = 0;
%% Read the specified File
fid = fopen(filename);
tline = fgetl(fid);
while ischar(tline)
    if (size(tline, 2) ~= 0 && tline(1) ~= '#')
        pattern = 'Constants';
        TF = contains(tline, pattern);
        if (TF)
            fprintf(1,'%s\n',tline);
            pos = strfind(tline, '=');
            exp = extractAfter(tline, pos(1));
            [C, ~] = strsplit(exp, {';'},'CollapseDelimiters',true);
            constants = sym(char(C(1)));
            flag_constants = 1;
        else
            [C, ~] = strsplit(tline,{'=', ';'},'CollapseDelimiters',true);
            exp = string(C(1));
            exp = strtrim(exp);
            switch exp
                case 'StateVariables'
                    fprintf(1,'%s\n',tline);
                    sv = transpose(sym(char(C(2))));
                    nx = size(sv, 1);
                case 'ControlVariables'
                    fprintf(1,'%s\n',tline);
                    cv = transpose(sym(char(C(2))));
                    nu = size(cv, 1);
                case 'ParameterVariables'
                    fprintf(1,'%s\n',tline);
                    pv = transpose(sym(char(C(2))));
                    nw = size(pv, 1);
                case 'InitialConstraints'
                    fprintf(1,'%s\n',tline);
                    ic = sym(char(C(2)));
                    nT = size(ic, 2);
                case 'TerminalConstraints'
                    fprintf(1,'%s\n',tline);
                    tc = sym(char(C(2)));
                    nF = size(tc, 2);
                case 'TerminalPenalty'
                    fprintf(1,'%s\n',tline);
                    tp = sym(char(C(2)));
                    ntp = size(tp, 2);
                case 'CostFunctional'
                    fprintf(1,'%s\n',tline);
                    cf = sym(char(C(2)));
                case 'DifferentialEquations'
                    fprintf(1,'%s\n',tline);
                    de = sym(char(C(2)));
                case 'InequalityConstraints'
                    fprintf(1,'%s\n',tline);
                    cvic = sym(char(C(2)));
                    nd = size(cvic, 2);
                case 'EqualityConstraints'
                    fprintf(1,'%s\n',tline);
                    ec = sym(char(C(2)));
                    nau = size(ec, 2);
                case 'StateVariableInequalityConstraints'
                    fprintf(1,'%s\n',tline);
                    svic = sym(char(C(2)));
                    ns = size(svic, 2);
                case 'Nodes'
                    fprintf(1,'%s\n',tline);
                    Nodes = str2double(char(C(2)));
                case 'Tolerance'
                    fprintf(1,'%s\n',tline);
                    tol = str2double(char(C(2)));
                case 'OutputFile'
                    fprintf(1,'%s\n',tline);
                    exp = erase(char(C(2)), '.');
                    exp = strtrim(exp);
                    output_File = strip(exp, 'both', '"');
                case 'InputFile'
                    fprintf(1,'%s\n',tline);
                    exp = erase(char(C(2)), '.');
                    exp = strtrim(exp);
                    input_File = strip(exp, 'both', '"');
                    flag_inputFile = 1;
                case 'StateEstimate'
                    fprintf(1,'%s\n',tline);
                    state_estimate = sym(char(C(2)));
                    flag_se = 1;
                case 'ControlEstimate'
                    fprintf(1,'%s\n',tline);
                    control_estimate = sym(char(C(2)));
                    flag_ce = 1;
                case 'ParameterEstimate'
                    fprintf(1,'%s\n',tline);
                    parameter_estimate = sym(char(C(2)));
                    flag_pe = 1;
                case 'InitialTime'
                    fprintf(1,'%s\n',tline);
                    t0 = sym(char(C(2)));
                case 'FinalTime'
                    fprintf(1,'%s\n',tline);
                    tf = sym(char(C(2)));
                case 'MaximumNodes'
                    fprintf(1,'%s\n',tline);
                    nodes_max = str2double(char(C(2)));
                case 'MaximumMeshRefinements'
                    fprintf(1,'%s\n',tline);
                    mesh_max = str2double(char(C(2)));
            end
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
if (nw == 0)
    pv = sym('pv', [nw, 1]);
end
if (nau == 0)
    ec = sym('ec', [1, nau]);
end
if (nd == 0)
    cvic = sym('cvic', [1, nd]);
end
if (ns == 0)
    svic = sym('svic', [1, ns]);
end
if (ntp == 0)
    tp = 0;
end
if (nd || ns)
    alpham = 1e-6;
end
%% Substitue the constants term
if (flag_constants)
    for i = 1 : size(constants, 2)
        ic = subs(ic, lhs(constants(i)), rhs(constants(i)));
        tc = subs(tc, lhs(constants(i)), rhs(constants(i)));
        cf = subs(cf, lhs(constants(i)), rhs(constants(i)));
        de = subs(de, lhs(constants(i)), rhs(constants(i)));
        t0 = subs(t0, lhs(constants(i)), rhs(constants(i)));
        tf = subs(tf, lhs(constants(i)), rhs(constants(i)));
        if (nau ~= 0)
            ec  = subs(ec, lhs(constants(i)), rhs(constants(i)));
        end
        if (nd ~= 0)
            cvic  = subs(cvic, lhs(constants(i)), rhs(constants(i)));
        end
        if (ns ~= 0)
            svic  = subs(svic, lhs(constants(i)), rhs(constants(i)));
        end
        if (ntp ~= 0)
            tp  = subs(tp, lhs(constants(i)), rhs(constants(i)));
        end
        if (flag_se)
            state_estimate = subs(state_estimate, lhs(constants(i)), rhs(constants(i)));
        end
        if (flag_ce)
            control_estimate = subs(control_estimate, lhs(constants(i)), rhs(constants(i)));
        end
        if (flag_pe)
            parameter_estimate = subs(parameter_estimate, lhs(constants(i)), rhs(constants(i)));
        end
    end
end

%% Generate the symbolic functions needed
ny = 2 * nx + nw;
nz = nu + 2 * nd + 2 * ns + nau;
np = nw + nT + nF;

% Define the variables
t = sym('t'); % time variable
alpha = sym('alpha'); % continuation parameter
lambda = sym('lambda', [nx 1]); % co-state variable
gama = sym('gama', [nw 1]); % co-parameter vector
mu = sym('mu', [nd 1]); % Lagrange multipliers a.w. the mixed control-state-parameter inequality constraints
e = sym('e', [ns 1]); % Lagrange multipliers a.w. the state variable inequality constraints
v = sym('v', [nd 1]); % non-negative slack variables a.w. the CVIC
sigma = sym('sigma', [ns 1]); % non-negative slack variables a.w. the SVIC
au = sym('au', [nau 1]); % Lagrange multipliers a.w. the auxiliary equality constraints
Ki = sym('Ki', [nT 1]); % Lagrange multipliers a.w. the initial constraints
Kf = sym('Kf', [nF 1]); % Lagrange multipliers a.w. the terminal constraints

y = [sv; lambda; gama];
z = [au; cv; mu; e; v; sigma];
p = [pv; Ki; Kf];

% Hamiltonian function
% H = 0.5 * alpha * transpose(cv) * cv + cf + de * lambda + ec * au + cvic * mu + svic * e;
H = cf + de * lambda + ec * au + cvic * mu + svic * e;
if (nd || ns)
    H = H + 0.5 * alpha * transpose(cv) * cv;
end
matlabFunction(H,'File','hamiltonian','Optimize', false, 'Vars', {[y; z], p, alpha});

%% DAEs of the system
h = sym('h', [ny 1]); % ODEs
g = sym('g', [nz 1]); % DAEs

% H / lamda
for i = 1 : nx
    h(i) = diff(H, lambda(i));
end
% H / x
for i = nx + 1 : 2 * nx
    h(i) = diff(-H, sv(i - nx));
end
% H / w
for i = 2 * nx + 1 : ny
    h(i) = diff(H, pv(i - 2 * nx));
end
% H / au
for i = 1 : nau
    g(i) = diff(H, au(i));
end
% H /u
for i = 1 : nu
    g(i + nau) = diff(H, cv(i));
end

%% Inequality constraint part
M = diag(mu);
N = diag(v);
E = diag(e);
S = diag(sigma);
e_d = ones(nd, 1);
e_s = ones(ns, 1);
g(nau + nu + 1 : nau + nu + nd) = transpose(cvic) + N * e_d - alpha * M * e_d;
g(nau + nu + nd + 1 : nau + nu + nd + ns) = transpose(svic) + S * e_s - alpha * E * e_s;

% Fischer-Burmeister formula
syms FBa FBb
psi(FBa, FBb) =  FBa + FBb - sqrt(FBa ^ 2 + FBb ^2 + 2 * alpha);
g(nau + nu + nd + ns + 1 : nau + nu + 2 * nd + ns) = psi(mu, v);
g(nau + nu + 2 * nd + ns + 1 : nau + nu + 2 * nd + 2 * ns) = psi(e, sigma);

matlabFunction(g,'File','DAE_g','Optimize', false, 'Vars', {y, z, p, alpha});
matlabFunction(h,'File','ODE_h','Optimize', false, 'Vars', {y, z, p, alpha});
%% Derivative of the DAEs
hy = sym('hy_%d_%d', [ny ny]);
hz = sym('hz_%d_%d', [ny nz]);
hp = sym('hp_%d_%d', [ny np]);
gy = sym('gy_%d_%d', [nz ny]);
gz = sym('gz_%d_%d', [nz nz]);
gp = sym('gp_%d_%d', [nz np]);

for i = 1 : ny
    for j = 1: ny
        hy(i, j) = diff(h(i), y(j));
    end
end
for i = 1 : ny
    for j = 1: nz
        hz(i, j) = diff(h(i), z(j));
    end
end
for i = 1 : ny
    for j = 1: np
        hp(i, j) = diff(h(i), p(j));
    end
end
for i = 1 : nz
    for j = 1: ny
        gy(i, j) = diff(g(i), y(j));
    end
end
for i = 1 : nz
    for j = 1: nz
        gz(i, j) = diff(g(i), z(j));
    end
end
for i = 1 : nz
    for j = 1: np
        gp(i, j) = diff(g(i), p(j));
    end
end
J1 = [hy hz
    gy gz];

matlabFunction(hy, hz, hp, gy, gz, gp,'File','difference_DAE','Optimize', false, 'Vars', {[y; z], p, alpha}, 'Sparse', true);
matlabFunction(J1,'File','jacobian_DAE','Optimize', false, 'Vars', {t, [y; z], p, alpha}, 'Sparse', true);

%% Define boundary constraints
sv0 = sym('sv%d_0', [nx 1]);
svN = sym('sv%d_N', [nx 1]);
lambda0 = sym('lambda%d_0', [nx 1]);
lambdaN = sym('lambda%d_N', [nx 1]);
gama0 = sym('gama%d_0', [nw 1]);
gamaN = sym('gama%d_N', [nw 1]);

r = sym('r', [ny + np 1]); % boundary constraints
D_Gamma_x = sym('D_Tao_x', [nx nT]);
D_Psi_x = sym('D_Psi_x', [nx nF]);
D_Phi_x = sym('D_fai_x', [nx, 1]);
D_Gamma_w = sym('D_Tao_w', [nw nT]);
D_Psi_w = sym('D_Psi_w', [nw nF]);
D_Phi_w = sym('D_fai_w', [nw, 1]);

y0 = [sv0; lambda0; gama0];
yM = [svN; lambdaN; gamaN];

% substitute the variable to distinguish the variables at the initial and
% terminal time
ic = subs(ic, sv, sv0);
tc = subs(tc, sv, svN);
tp = subs(tp, sv, svN);

for i = 1 : nT
    r(i) = ic(i);
end
for i = 1 : nx
    for j = 1 : nT
        D_Gamma_x(i, j) = diff(ic(j), sv0(i));
    end
end
r(1 + nT : nx + nT) = lambda0 + D_Gamma_x * Ki;
for i = 1 : nw
    for j = 1 : nT
        D_Gamma_w(i, j) = diff(ic(j), pv(i));
    end
end
r(1 + nT + nx : nw + nT + nx) = gama0 - D_Gamma_w * Ki;
for i = 1 : nF
    r(i + nT + nx + nw) = tc(i);
end
for i = 1 : nx
    D_Phi_x(i) = diff(tp, svN(i));
end
for i = 1 : nx
    for j = 1 : nF
        D_Psi_x(i, j) = diff(tc(j), svN(i));
    end
end
r(1 + nT + nx + nw + nF : nx + nT + nx + nw + nF) = lambdaN - D_Phi_x - D_Psi_x * Kf;
for i = 1 : nw
    D_Phi_w(i) = diff(tp, pv(i));
end
for i = 1 : nw
    for j = 1 : nF
        D_Psi_w(i, j) = diff(tc(j), pv(i));
    end
end
r(1 + nT + nx + nw + nF + nx : nw + nT + nx + nw + nF + nx) = gamaN + D_Phi_w + D_Psi_w * Kf;
matlabFunction(r,'File','boundary_constraint','Optimize', false, 'Vars', {y0, yM, p});

%% Derivative of bounday constraints
r_x0 = sym('r_x0', [ny + np ny]);
r_xM = sym('r_xM', [ny + np ny]);
r_p = sym('r_p', [ny + np np]);
% r w.r.t y0
for i = 1 : ny + np
    for j = 1 : ny
        r_x0(i, j) = diff(r(i), y0(j));
    end
end
% r w.r.t yM
for i = 1 : ny + np
    for j = 1 : ny
        r_xM(i, j) = diff(r(i), yM(j));
    end
end
% r w.r.t p
for i = 1 : ny + np
    for j = 1: np
        r_p(i, j) = diff(r(i), p(j));
    end
end
matlabFunction(r_x0,r_xM,r_p,'File','difference_BC','Optimize', false, 'Vars', {y0, yM, p});

%% Set up state and control estimate
t0 = double(t0);
tf = double(tf);
t_span = linspace(t0, tf, Nodes)';
y0 = ones(Nodes, ny);
z0 = ones(Nodes, nz);

if (flag_se || flag_ce)
    for i = 1 : Nodes
        if (flag_se)
            for j = 1 : size(state_estimate, 2)
                y0(i, j) = subs(state_estimate(j), t, t_span(i));
            end
        end
        if (flag_ce)
            for j = 1 : size(control_estimate, 2)
                z0(i, j) = subs(control_estimate(j), t, t_span(i));
            end
        end
    end
elseif (flag_inputFile)
    input_File_y = [input_File, '_y.txt'];
    input_File_z = [input_File, '_z.txt'];
    input_y = load(input_File_y);
    input_z = load(input_File_z);
    size_y = size(input_y, 2);
    Nodes = size(input_y, 1);
    y0 = ones(Nodes, ny);
    z0 = ones(Nodes, nz);
    t_span = input_y(1 : Nodes, 1);
    y0(1 : Nodes, 1 : size_y - 1) = input_y(1 : Nodes, 2 : end);
    z0(1 : Nodes, 1 : nu) = input_z(1 : Nodes, 1 : nu);
    z0(Nodes, 1 : nu) = input_z(end, 1 : nu);
end
% y0(1 : Nodes, 1 + nx : ny) = 10 * ones(Nodes, (ny - nx));
% z0(1 : Nodes, 1 + nu : nz) = 10 * ones(Nodes, (nz - nu));
para = ones(np, 1);
if (flag_pe)
    for j = 1 : size(parameter_estimate, 2)
        para(j) = double(parameter_estimate(j));
    end
end
end