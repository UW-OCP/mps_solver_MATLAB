function mps_BVP(ODE_h, DAE_g, BC_r, D_hg, D_r, J_hg, y0, z0, p, ti, tM, output, alpham, tolerance, nodes_max, mesh_max, nx, nu, nw, nd, ns)

warning off
t_all = tic;

global ny nz np M tspan t0 tf tol coefficients

if (size(y0, 1) ~= size(z0, 1))
    error('Error. \nInput size is not right, input size of y (%d) is not equal to input size of z (%d).',size(y0, 1),size(z0, 1));
end

ny = size(y0, 2);
nz = size(z0, 2);
np = size(p, 1);
M = size(y0, 1);
t0 = ti;
tf = tM;
tspan = linspace((t0),(tf),M);

tol = tolerance;
MAX_ITER = 500;
Maxerr = 1e12;
alpha0search = 20;

alpha = 1;
beta = 0.9;
flag = 0;
flag_coefficients = 1;

time_residual = 0;
time_construct = 0;
time_BABD = 0;
number_of_times_residual = 0;
number_of_times_construct = 0;
number_of_times_BABD = 0;

% file to print the run-time information
filename_info = [output, '_info.txt'];
fileID_info = fopen(filename_info,'w');
fprintf(fileID_info, '%s\n', 'Smoothed Fischer-Burmeister noninterior.');
fprintf(fileID_info, '%s%d\n', 'method: Multiple shooting. Flag coefficients: ', flag_coefficients);
%% Initial Guess
s0 = matrixtovec(y0, z0);
%% Get Integrator Coefficients
coefficients = getCoefficients(flag_coefficients);
%% Continuation Method
for alphacal = 1:MAX_ITER
    %% Multiple Shooting method
    mesh_time = 0;
    for caltime=1:MAX_ITER
        %% Computer F(s0)
        t_construct = tic;
        sol = construct(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
            @(s, p, alpha) D_hg(s, p, alpha), @(t, s, p, alpha) J_hg(t, s, p, alpha),...
            s0, p, alpha, tspan, coefficients);
        time_construct_it = toc(t_construct);
        time_construct = time_construct + time_construct_it;
        number_of_times_construct = number_of_times_construct + 1;
        
        t_residual = tic;
        [F_s0, sol, maxlte] = F_s_residual(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
            @(t, s, p, alpha) J_hg(t, s, p, alpha), @(y0, yM, p) BC_r(y0, yM, p),...
            sol, s0, p, tspan, alpha, 'yes');
%         [F_s0, sol, maxlte] = F_s_spmd_residual(sol, s0, p, tspan, alpha, 'yes');
        time_residual_it = toc(t_residual);
        time_residual = time_residual + time_residual_it;
        number_of_times_residual = number_of_times_residual + 1;
        
        if norm(F_s0, Inf) < tol
            disp(['alpha = ', num2str(alpha), '; nodes = ', num2str(M), '; LTE = ', num2str(maxlte), '.']);
            fprintf(fileID_info, '%s%d%s%d%s%d%s\n', 'alpha = ', alpha, '; nodes = ', M, '; LTE = ', maxlte, '.');
            break
        end
        if (caltime == MAX_ITER || norm(F_s0, Inf) > Maxerr)
            disp('Solution is not converging.')
            fprintf(fileID_info, '%s\n', 'Solution is not converging.');
            flag = 1;
            break
        end
        %% Generation of the Jacobian Matrix
        sol = DF(@(y0, yM, p) D_r(y0, yM, p), sol, p);
        t_BABD = tic;
%         delta_s = qr_spmd(sol);
        %% QR decomposition
        sol = sequential_qr(sol);
        %% Compute delta_s
        delta_s = backwardsubstitution(sol);
        time_BABD_it = toc(t_BABD);
        time_BABD = time_BABD + time_BABD_it;
        number_of_times_BABD = number_of_times_BABD + 1;
        if norm(delta_s(1 : end - np), Inf) < tol
            disp(['alpha = ', num2str(alpha), '; nodes = ', num2str(M), '; LTE = ', num2str(maxlte), '.']);
            fprintf(fileID_info, '%s%d%s%d%s\n', 'alpha = ', alpha, '; nodes = ', M, '.');
            break
        end
        
        %% Find an alpha linesearch
        if (isnan(norm(delta_s)))
            [s0, tspan] = mesh_refinement(sol, s0, tspan);
            mesh_time = mesh_time + 1;
            disp(['Wrong Direction! Need to remesh!, Number of nodes = ', num2str(M), '.']);
            fprintf(fileID_info, '%s%d%s\n', 'Wrong Direction! Need to remesh!, Number of nodes = ', M, '.');
        else
            alpha0 = 1;
            for i=1 : alpha0search
                s_New = s0 + alpha0 * delta_s(1:end - np);
                p_New = p + alpha0 * delta_s(1 + end - np : end);
                
                t_residual = tic;
                [F_s0_new, sol, ~] = F_s_residual(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
                    @(t, s, p, alpha) J_hg(t, s, p, alpha), @(y0, yM, p) BC_r(y0, yM, p),...
                    sol, s_New, p_New, tspan, alpha, 'no');
%                 [F_s0_new, sol, ~] = F_s_spmd_residual(sol, s_New, p_New, tspan, alpha, 'no');
                time_residual_it = toc(t_residual);
                time_residual = time_residual + time_residual_it;
                number_of_times_residual = number_of_times_residual + 1;
                
                if norm(F_s0_new, Inf) < norm(F_s0, Inf) || norm(F_s0_new, Inf) < tol
                    s0 = s_New;
                    p = p_New;
                    break
                end
                alpha0 = alpha0/2;
            end
            if i == alpha0search
                disp(['The solution does not converge at the alpha line-search, Remesh!, alpha = ', num2str(alpha), '.']);
                fprintf(fileID_info, '%s%d%s%d%s\n', 'The solution does not converge at the alpha line-search, Remesh!, alpha = ', alpha, '; Nodes = ', M, '.');
                [s0, tspan] = mesh_refinement(sol, s0, tspan);
                mesh_time = mesh_time + 1;
                disp(['Remeshed!, Number of nodes = ', num2str(M), '.']);
                fprintf(fileID_info, '%s%d%s\n', 'Remeshed!, Number of nodes = ', M, '.');
            end
        end
        if M >= nodes_max || M <= 20
            disp(['Number of nodes is too many or too few! Nodes = ', num2str(M), '.']);
            fprintf(fileID_info, '%s%d%s\n', 'Number of nodes is too many or too few! Nodes = ', M, '.');
            flag = 1;
            break;
        elseif mesh_time > mesh_max
            disp('Failed to find the right solutino after too many meshes!');
            fprintf(fileID_info, '%s\n', 'Failed to find the right solutino after too many meshes!');
            flag = 1;
            break;
        end
    end
    if flag == 1
        disp(['The solution does not converge, Stop!, alpha = ', num2str(alpha), '.']);
        fprintf(fileID_info, '%s%d%s\n', 'The solution does not converge, Stop!, alpha = ', alpha, '.');
        break
    end
    while maxlte > 1
        [s0, tspan] = mesh_refinement(sol, s0, tspan);
        t_residual = tic;
        [~, sol, maxlte] = F_s_residual(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
            @(t, s, p, alpha) J_hg(t, s, p, alpha), @(y0, yM, p) BC_r(y0, yM, p),...
            sol, s0, p, tspan, alpha, 'yes');
%         [~, sol, maxlte] = F_s_spmd_residual(sol, s0, p, tspan, alpha, 'yes');
        time_residual_it = toc(t_residual);
        time_residual = time_residual + time_residual_it;
        number_of_times_residual = number_of_times_residual + 1;
        disp(['LTE = ', num2str(maxlte), '. Remesh! Number of Nodes = ', num2str(M), '.']);
        fprintf(fileID_info, '%s%d%s%d%s\n', 'LTE = ', maxlte, '. Remesh! Number of Nodes = ', M, '.');
    end
    if alpha <= alpham
        disp(['Get the final solution!, alpha = ', num2str(alpha), '.'])
        fprintf(fileID_info, '%s%d%s\n', 'Get the final solution! alpha = ', alpha);
        break
    end
    alpha = beta*alpha;
    if alpha < tol
        tol = 0.9*alpha;
    end
end
time_all = toc(t_all);
fclose(fileID_info);
%% Save Final Result
[sol_y, sol_z] = vectomatrix(s0);
filename_y = [output, '_y.txt'];
fileID_y = fopen(filename_y,'w');
spec = '%8.6f';
for i = 1 : ny
    spec = [spec, ' %12.8f'];
end
spec = [spec, '\n'];
sol_y_print = [tspan' sol_y]';
fprintf(fileID_y, spec, sol_y_print);
fclose(fileID_y);
filename_z = [output, '_z.txt'];
fileID_z = fopen(filename_z,'w');
spec = [];
for i = 1 : nz
    spec = [spec, ' %12.8f'];
end
spec = [spec, '\n'];
fprintf(fileID_z, spec, (sol_z)');
fclose(fileID_z);
filename_p = [output, '_p.txt'];
fileID1 = fopen(filename_p,'w');
fprintf(fileID1,'%12.8f\n',p);
fclose(fileID1);
filename_alpha = [output, '_alpha.txt'];
fileID_alpha = fopen(filename_alpha,'w');
fprintf(fileID_alpha,'%12.8f\n',alpha);
fclose(fileID_alpha);
filename_time = [output, '_time.txt'];
fileID_time = fopen(filename_time, 'w');
fprintf(fileID_time,'%12.8f\n',time_all);
fclose(fileID_time);
filename_time_construct = [output, '_time_construct.txt'];
fileID_time_construct = fopen(filename_time_construct, 'w');
fprintf(fileID_time_construct,'%12.8f\n',time_construct);
fclose(fileID_time_construct);
filename_time_residual = [output, '_time_residual.txt'];
fileID_time_residual = fopen(filename_time_residual, 'w');
fprintf(fileID_time_residual,'%12.8f\n',time_residual);
fclose(fileID_time_residual);
filename_time_BABD = [output, '_time_BABD.txt'];
fileID_time_BABD = fopen(filename_time_BABD, 'w');
fprintf(fileID_time_BABD,'%12.8f\n',time_BABD);
fclose(fileID_time_BABD);

filename_times_construct = [output, '_times_construct.txt'];
fileID_times_construct = fopen(filename_times_construct, 'w');
fprintf(fileID_times_construct,'%d\n',number_of_times_construct);
fclose(fileID_times_construct);
filename_times_residual = [output, '_times_residual.txt'];
fileID_times_residual = fopen(filename_times_residual, 'w');
fprintf(fileID_times_residual,'%d\n',number_of_times_residual);
fclose(fileID_times_residual);
filename_times_BABD = [output, '_times_BABD.txt'];
fileID_times_BABD = fopen(filename_times_BABD, 'w');
fprintf(fileID_times_BABD,'%d\n',number_of_times_BABD);
fclose(fileID_times_BABD);

%% Plot the results
t_span = tspan;
for i = 1 : nx
    figure();
    plot(t_span, sol_y(:, i), 'DisplayName', strcat('x_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('State');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : nx
    figure();
    plot(t_span, sol_y(:, i + nx), 'DisplayName', strcat('\lambda_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Costate');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : nw
    figure();
    plot(t_span, sol_y(:, i + 2 * nx), 'DisplayName', strcat('\gamma_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Co-Parameter-State');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : nu
    figure();
    plot(t_span(1 : end), sol_z(1 : end, i), 'DisplayName', strcat('u_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Control');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : nd
    figure();
    plot(t_span(1 : end), sol_z(1 : end, nu + i), 'DisplayName', strcat('\mu_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Lagrange Parameter w.r.t CVIC');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : ns
    figure();
    plot(t_span(1 : end), sol_z(1 : end, nu + nd + i), 'DisplayName', strcat('\sigma_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Lagrange Parameter w.r.t SVIC');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : nd
    figure();
    plot(t_span(1 : end), sol_z(1 : end, nu + nd + ns + i), 'DisplayName', strcat('v_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Lagrange Parameter w.r.t CVIC');
    grid on;
    legend('Location', 'best');
    legend show;
end
for i = 1 : ns
    figure();
    plot(t_span(1 : end), sol_z(1 : end, nu + 2 * nd + ns + i), 'DisplayName', strcat('\xi_{', num2str(i), '}'), 'LineWidth', 1);
    xlabel('Time');
    ylabel('Lagrange Parameter w.r.t SVIC');
    grid on;
    legend('Location', 'best');
    legend show;
end
end

%% Coefficients of the RK
function coefficients = getCoefficients(flag)

if flag == 1
    %schneider
    s = 8;
    stages = s;
    order = 5;
    alpha = zeros(s,s);
    beta = zeros(s,s);
    gamma = zeros(s,s);
    b = zeros(s, 1);
    e = zeros(s, 1);
    
    alpha(2, 1) = 0.5;
    
    alpha(3, 1) = 0.25;
    alpha(3, 2) = 0.25;
    
    alpha(4, 1) = 0.25;
    alpha(4, 2) = 0.25;
    
    alpha(5, 1) = 9.0/32.0;
    alpha(5, 2) = -21.0/32;
    alpha(5, 3) = -117.0/184.0;
    alpha(5, 4) = 81.0/46.0;
    
    alpha(6, 1) = 97.0/144.0;
    alpha(6, 2) = -751.0/1536.0;
    alpha(6, 3) = 26107.0/11776.0;
    alpha(6, 4) = -3249.0/1472.0;
    alpha(6, 5) = 29.0/36.0;
    
    alpha(7, 1) = 5.0/27.0;
    alpha(7, 2) = -65.0/128.0;
    alpha(7, 3) = -3873.0/2944.0;
    alpha(7, 4) = 947.0/368.0;
    alpha(7, 5) = -5.0/27.0;
    alpha(7, 6) = 0.25;
    
    alpha(8, 1) = 71.0/270.0;
    alpha(8, 4) = 2.0/5.0;
    alpha(8, 5) = 32.0/135.0;
    alpha(8, 6) = -3.0/20.0;
    alpha(8, 7) = 0.25;
    
    for i = 1 : s
        beta(i, i) = 0.25;
    end
    
    beta(3, 1) = -1.0/12.0;
    beta(3, 2) = -1.0/12.0;
    
    beta(4, 1) = 11.0/114.0;
    beta(4, 2) = -13.0/64.0;
    beta(4, 3) = 23.0/192.0;
    
    beta(5, 1) = 0.75;
    beta(5, 2) = -819.0/1024.0;
    beta(5, 3) = 11061.0/23552.0;
    beta(5, 4) = 243.0/736.0;
    
    beta(6, 1) = 5.0/27.0;
    beta(6, 2) = -65.0/128.0;
    beta(6, 3) = -3873.0/2944.0;
    beta(6, 4) = 947.0/368.0;
    beta(6, 5) = -5.0/27.0;
    
    beta(7, 1) = 71.0/270.0;
    beta(7, 4) = 2.0/5.0;
    beta(7, 5) = 32.0/135.0;
    beta(7, 6) = -3.0/20.0;
    
    beta(8, 1) = 71.0/270.0;
    beta(8, 4) = 2.0/5.0;
    beta(8, 5) = 32.0/135.0;
    beta(8, 6) = -8.0/15.0;
    beta(8, 7) = 23.0/60.0;
    
    for i = 1 : s
        b(i) = alpha(s, i);
        for j = 1 : i
            gamma(i, j) = beta(i, j) - alpha(i, j);
        end
    end
    
    e(6) = -8.0/15.0 + 3.0/20.0;
    e(7) = 23.0/60.0 - 0.25;
    e(8) = 0.25;
elseif flag == 2
    % rodas
    gamma_ = 0.25;
    s = 6;
    stages = s;
    order = 4;
    alpha = zeros(s,s);
    gamma = zeros(s,s);
    b = zeros(s);
    e = zeros(s);
    for i = 1 : s
        gamma(i, i) = gamma_;
    end
    a21 = 0.386;
    a31 = 0.146074707525418;
    a32 = 0.063925292474582;
    a41 = -0.330811503667722;
    a42 = 0.711151025168282;
    a43 = 0.24966047849944;
    a51 = -4.552557186318003;
    a52 = 1.710181363241322;
    a53 = 4.014347332103150;
    a54 = -0.171971509026469;
    a61 = 2.428633765466978;
    a62 = -0.382748733764781;
    a63 = -1.855720330929574;
    a64 = 0.559835299227375;
    a65 = 0.25;
    g21 = -0.3543;
    g31 = -0.133602505268175;
    g32 = -0.012897494731825;
    g41 = 1.526849173006459;
    g42 = -0.533656288750454;
    g43 = -1.279392884256;
    g51 = 6.981190951784981;
    g52 = -2.092930097006103;
    g53 = -5.870067663032724;
    g54 = 0.731806808253845;
    g61 = -2.080189494180926;
    g62 = 0.59576235567668;
    g63 = 1.701617798267255;
    g64 = -0.088514519835879;
    g65 = -0.378676139927128;
    b1 = 0.348444271286054;
    b2 = 0.213013621911897;
    b3 = -0.154102532662319;
    b4 = 0.471320779391497;
    b5 = -0.128676139927129;
    b6 = 0.25;
    be1 = -4.552557186318003;
    be2 = 1.710181363241322;
    be3 = 4.014347332103150;
    be4 = -0.171971509026469;
    be5 = 0;
    be6 = 0;
    alpha(2, 1) = a21;
    alpha(3, 1) = a31;
    alpha(3, 2) = a32;
    alpha(4, 1) = a41;
    alpha(4, 2) = a42;
    alpha(4, 3) = a43;
    alpha(5, 1) = a51;
    alpha(5, 2) = a52;
    alpha(5, 3) = a53;
    alpha(5, 4) = a54;
    alpha(6, 1) = a61;
    alpha(6, 2) = a62;
    alpha(6, 3) = a63;
    alpha(6, 4) = a64;
    alpha(6, 5) = a65;
    gamma(2, 1) = g21;
    gamma(3, 1) = g31;
    gamma(3, 2) = g32;
    gamma(4, 1) = g41;
    gamma(4, 2) = g42;
    gamma(4, 3) = g43;
    gamma(5, 1) = g51;
    gamma(5, 2) = g52;
    gamma(5, 3) = g53;
    gamma(5, 4) = g54;
    gamma(6, 1) = g61;
    gamma(6, 2) = g62;
    gamma(6, 3) = g63;
    gamma(6, 4) = g64;
    gamma(6, 5) = g65;
    b(1) = b1;
    b(2) = b2;
    b(3) = b3;
    b(4) = b4;
    b(5) = b5;
    b(6) = b6;
    e(1) = b1 - be1;
    e(2) = b2 - be2;
    e(3) = b3 - be3;
    e(4) = b4 - be4;
    e(5) = b5 - be5;
    e(6) = b6 - be6;
elseif flag == 3
    % rodasp
    gamma_ = 0.25;
    s = 6;
    stages = s;
    order = 4;
    alpha = zeros(s,s);
    gamma = zeros(s,s);
    b = zeros(s);
    e = zeros(s);
    for i = 1 : s
        gamma(i, i) = gamma_;
    end
    alpha(2, 1) = 7.500000e-01;
    alpha(3, 1) = 8.612040e-02;
    alpha(3, 2) = 1.238796e-01;
    alpha(4, 1) = 7.749345e-01;
    alpha(4, 2) = 1.492652e-01;
    alpha(4, 3) = -2.941997e-01;
    alpha(5, 1) = 5.308747e+00;
    alpha(5, 2) = 1.330892e+00;
    alpha(5, 3) = -5.374138e+00;
    alpha(5, 4) = -2.655010e-01;
    alpha(6, 1) = -1.764438e+00;
    alpha(6, 2) = -4.747566e-01;
    alpha(6, 3) = 2.369692e+00;
    alpha(6, 4) = 6.195024e-01;
    alpha(6, 5) = 2.500000e-01;
    gamma(1, 1) = 2.500000e-01;
    gamma(2, 1) = -7.500000e-01;
    gamma(2, 2) = 2.500000e-01;
    gamma(3, 1) = -1.355124e-01;
    gamma(3, 2) = -1.379916e-01;
    gamma(3, 3) = 2.500000e-01;
    gamma(4, 1) = -1.256984e+00;
    gamma(4, 2) = -2.501447e-01;
    gamma(4, 3) = 1.220929e+00;
    gamma(4, 4) = 2.500000e-01;
    gamma(5, 1) = -7.073184e+00;
    gamma(5, 2) = -1.805649e+00;
    gamma(5, 3) = 7.743830e+00;
    gamma(5, 4) = 8.850034e-01;
    gamma(5, 5) = 2.500000e-01;
    gamma(6, 1) = 1.684069e+00;
    gamma(6, 2) = 4.182659e-01;
    gamma(6, 3) = -1.881406e+00;
    gamma(6, 4) = -1.137861e-01;
    gamma(6, 5) = -3.571429e-01;
    gamma(6, 6) = 2.500000e-01;
    b(1) = -8.036837e-02;
    b(2) = -5.649061e-02;
    b(3) = 4.882856e-01;
    b(4) = 5.057162e-01;
    b(5) = -1.071429e-01;
    b(6) = 2.500000e-01;
    e(1) = 1.684069e+00;
    e(2) = 4.182659e-01;
    e(3) = -1.881406e+00;
    e(4) = -1.137861e-01;
    e(5) = -3.571429e-01;
    e(6) = 2.500000e-01;
elseif flag == 4
    % wb34
    gamma_ = 5.728160624821350e-01;
    s = 6;
    stages = s;
    order = 4;
    alpha = zeros(s,s);
    gamma = zeros(s,s);
    b = zeros(s, 1);
    e = zeros(s);
    for i = 1 : s
        gamma(i, i) = gamma_;
    end
    a21 = 5.200000000000000e-01;
    g21 = -5.200000000000000e-01;
    a31 = 2.851168665349716e-01;
    g31 = -1.034772479328808e+00;
    a32 = 6.248831334650284e-01;
    g32 = 6.501423878169246e-01;
    a41 = 1.046681454850720e+00;
    g41 = 2.625385974420247e-01;
    a42 = -1.127221164631929e+00;
    g42 = 2.922670258511625e-01;
    a43 = 3.910371962111624e-01;
    g43 = -9.114397095544884e-01;
    a51 = 8.451547656533995e-02;
    g51 = 1.574388804512719e-01;
    a52 = 1.140000000000000e+00;
    g52 = 6.277349506307095e-02;
    a53 = -6.668002390497316e-02;
    g53 = -5.710378229055593e-01;
    a54 = -1.578354526603668e-01;
    g54 = -2.219906150909184e-01;
    a61 = 2.419543570166118e-01;
    g61 = 0.000000000000000e-00;
    a62 = 1.202773495063071e+00;
    g62 = 0.000000000000000e-00;
    a63 = -6.377178468105325e-01;
    g63 = 0.000000000000000e-00;
    a64 = -3.798260677512852e-01;
    g64 = 0.000000000000000e-00;
    a65 = 5.728160624821350e-01;
    g65 = -5.728160624821350e-01;
    b1 = 2.419543570166118e-01;
    b1e = 2.419543570166118e-01;
    b2 = 1.202773495063071e+00;
    b2e = 1.202773495063071e+00;
    b3 = -6.377178468105325e-01;
    b3e = -6.377178468105325e-01;
    b4 = -3.798260677512852e-01;
    b4e = -3.798260677512852e-01;
    b5 = 0.000000000000000e-00;
    b5e = 5.728160624821350e-01;
    b6 = 5.728160624821350e-01;
    b6e = 0.000000000000000e-00;
    alpha(2, 1) = a21;
    alpha(3, 1) = a31;
    alpha(3, 2) = a32;
    alpha(4, 1) = a41;
    alpha(4, 2) = a42;
    alpha(4, 3) = a43;
    alpha(5, 1) = a51;
    alpha(5, 2) = a52;
    alpha(5, 3) = a53;
    alpha(5, 4) = a54;
    alpha(6, 1) = a61;
    alpha(6, 2) = a62;
    alpha(6, 3) = a63;
    alpha(6, 4) = a64;
    alpha(6, 5) = a65;
    gamma(2, 1) = g21;
    gamma(3, 1) = g31;
    gamma(3, 2) = g32;
    gamma(4, 1) = g41;
    gamma(4, 2) = g42;
    gamma(4, 3) = g43;
    gamma(5, 1) = g51;
    gamma(5, 2) = g52;
    gamma(5, 3) = g53;
    gamma(5, 4) = g54;
    gamma(6, 1) = g61;
    gamma(6, 2) = g62;
    gamma(6, 3) = g63;
    gamma(6, 4) = g64;
    gamma(6, 5) = g65;
    b(1) = b1;
    b(2) = b2;
    b(3) = b3;
    b(4) = b4;
    b(5) = b5;
    b(6) = b6;
    e(1) = b1 - b1e;
    e(2) = b2 - b2e;
    e(3) = b3 - b3e;
    e(4) = b4 - b4e;
    e(5) = b5 - b5e;
    e(6) = b6 - b6e;
end
coefficients.stages = stages;
coefficients.alpha = alpha;
coefficients.gamma = gamma;
coefficients.b = b;
coefficients.e = e;
end

%% Construct the initial structure
function sol = construct(ODE_h, DAE_g, D_hg, J_hg,s, p, alpha, tspan, coefficients)

global M ny nz np
sol(M)=struct('y',[],'z',[],'h_y',[],'h_z',[],'h_p',[],'g_y',[],'g_z',[],'g_p',[],'y_y',[],'y_z',[],'y_p',[],...
    'A',[], 'C',[], 'H',[], 'b', [], 'B', [], 'index', [],...
    'R', [], 'E', [], 'J', [], 'G', [], 'C_til', [], 'A_til', [], 'H_til', [], 'G_til', [], 'b_til', [], 'delta_s', [], 'delta_p', []);
% parfor i = 1:M
for i = 1:M
    sol(i).y = s(1 + (i - 1) * (nz + ny) : ny + (i - 1) * (nz + ny));
    sol(i).z = s(1 + ny + (i - 1) * (nz + ny) : i * (nz + ny));
end
% parfor i = 1 : M - 1
for i = 1 : M - 1
    [h_y, h_z, h_p, g_y, g_z, g_p] = D_hg([sol(i).y; sol(i).z], p, alpha);
    sol(i).h_y = h_y;
    sol(i).h_z = h_z;
    sol(i).h_p = h_p;
    sol(i).g_y = g_y;
    sol(i).g_z = g_z;
    sol(i).g_p = g_p;
    t_span = [tspan(i) tspan(i + 1)];
    [~, X_next] = row_sensitivity_step(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
        @(s, p, alpha) D_hg(s, p, alpha), @(t, s, p, alpha) J_hg(t, s, p, alpha),...
        t_span, [sol(i).y; sol(i).z], p, alpha, ny, nz, np, coefficients);
    sol(i).y_y = X_next(1 : ny, 1 : ny); 
    sol(i).y_z = X_next(1 : ny, 1 + ny : nz + ny);
    sol(i).y_p = X_next(1 : ny, 1 + ny + nz : nz + ny + np);
end
[h_y, h_z, h_p, g_y, g_z, g_p] = D_hg([sol(M).y; sol(M).z], p, alpha);
sol(M).h_y = h_y;
sol(M).h_z = h_z;
sol(M).h_p = h_p;
sol(M).g_y = g_y;
sol(M).g_z = g_z;
sol(M).g_p = g_p;
end

%% Calculate the residual of the problem
function [F_s, sol, maxlte] = F_s_residual(ODE_h, DAE_g, J_hg, BC_r, sol, s, p, tspan, alpha, flag_update)

global M ny nz np coefficients tol

F_s = zeros((ny + nz) * M + np, 1);
residuals(M) = struct('b', []);
lte_all = zeros(M - 1, 1);
%% Integration of the DAE
parfor i = 1 : M - 1
    y = s(1 + (i - 1) * (nz + ny) : ny + (i - 1) * (nz + ny));
    z = s(1 + ny + (i - 1) * (nz + ny) : i * (nz + ny));
    t_span = [tspan(i); (tspan(i) + tspan(i + 1))/2; tspan(i + 1)];
    G = DAE_g(y, z, p, alpha);
    [x_next, lte] = row_step(@(y, z, p, alpha) ODE_h(y, z, p, alpha), @(y, z, p, alpha) DAE_g(y, z, p, alpha),...
        @(t, s, p, alpha) J_hg(t, s, p, alpha),...
        t_span, [y; z], p, alpha, ny, nz, np, coefficients, tol);
    y_next = x_next(1 : ny);
    residuals(i).b = -[G; y_next - s(1 + i * (nz + ny) : ny + i * (nz + ny))];
    lte_all(i) = lte;
end
maxlte = norm(lte_all, inf);
y_0 = s(1 : ny);
y_M = s(1 + (M - 1) * (nz + ny) : ny + (M - 1) * (nz + ny));
z = s(1 + ny + (M - 1) * (nz + ny) : M * (nz + ny));
residuals(M).b = -[DAE_g(y_M, z, p, alpha); BC_r(y_0, y_M, p)];

% update the residuals and the ltes
if (strcmp(flag_update, 'yes'))
    for i = 1 : M - 1
        F_s(1 + (i - 1) * (nz + ny) : i * (nz + ny)) = -residuals(i).b;
        sol(i).b = residuals(i).b;
        sol(i).e = lte_all(i);
    end
    F_s(1 + (M - 1) * (ny + nz) : end) = -residuals(M).b;
    sol(M).b = residuals(M).b;
else
    for i = 1 : M - 1
        F_s(1 + (i - 1) * (nz + ny) : i * (nz + ny)) = -residuals(i).b;
    end
    F_s(1 + (M - 1) * (ny + nz) : end) = -residuals(M).b;
end
end

%% Construct the Jacobian of the problem
function sol = DF(D_r, sol, p)
% generation of the Jacobian matrix

global M ny nz np

% for the B0&BM part
[r_s0, r_sM, r_p] = D_r(sol(1).y, sol(M).y, p);

B_1 = zeros(nz + ny + np, ny + nz);
B_1(1 + nz : nz + ny + np, 1 : ny) = r_s0;
sol(1).B = B_1;

B_N = zeros(nz + ny + np, ny + nz);
B_N(1 : nz, 1 : ny) = sol(M).g_y;
B_N(1 : nz, 1 + ny : ny + nz) = sol(M).g_z;
B_N(1 + nz : nz + ny + np, 1 : ny) = r_sM;
sol(M).B = B_N;

H_N = zeros(nz + ny + np, np);
H_N(1 : nz, 1 : np) = sol(M).g_p;
H_N(1 + nz : nz + ny + np, 1 : np) = r_p;
sol(M).H = H_N;

% for the A&C&H parts
parfor i = 1 : M - 1
    A = [sol(i).g_y sol(i).g_z
        sol(i).y_y sol(i).y_z];
    H = [sol(i).g_p
        sol(i).y_p];
    C_D = -ones(ny,1);
    C_i = 1 + nz : ny + nz;
    C_j = 1 : ny;
    C = sparse(C_i,C_j,C_D,(ny+nz),(ny+nz));
    sol(i).A = A;
    sol(i).C = C;
    sol(i).H = H;
end
end

%% Parallel qr decomposition with spmd
function delta_s = qr_spmd(sol)
global ny nz np M

p = gcp;
num_workers = p.NumWorkers;
% num_workers = 4;

ns = ny + nz;

% make a copy of the struct from index 1 to index M - 1
% also index the struct
sol_reduced(M - 1) = struct('A', [], 'C', [],'H', [], 'b', [], 'index', []);
for i = 1 : M - 1
    sol_reduced(i).A = sol(i).A;
    sol_reduced(i).C = sol(i).C;
    sol_reduced(i).H = sol(i).H;
    sol_reduced(i).b = sol(i).b;
    sol_reduced(i).index = i;
end

sol_reduced_dis = distributed(sol_reduced);
spmd
   sol_local = getLocalPart(sol_reduced_dis); % Unique value on each worker 
   size_sol = size(sol_local, 2);
%    r_start = sol_local(1).index;
%    r_end = sol_local(size_sol).index;
   sol_update(size_sol) = struct('R', [], 'E', [], 'J', [], 'G', [], 'C_til', [], 'A_til', [], 'H_til', [], 'G_til', [], 'b_til', [], 'delta_s', []);
   sol_update(1).C_til = sol_local(1).C;
   sol_update(1).G_til = sol_local(1).A;
   sol_update(1).H_til = sol_local(1).H;
   sol_update(1).b_til = sol_local(1).b;
   for i = 1 : size_sol - 1
       [Q, R] = qr([sol_update(i).C_til; sol_local(i + 1).A]);
       sol_update(i).R = R(1 : ns, :);
       EC = Q' * [zeros(ns); sol_local(i + 1).C];
       sol_update(i).E = EC(1 : ns, 1 : ns);
       sol_update(i + 1).C_til = EC(1 + ns : 2 * ns, 1 : ns);
       GG = Q' * [sol_update(i).G_til; zeros(ns)];
       sol_update(i).G = GG(1 : ns, 1 : ns);
       sol_update(i + 1).G_til = GG(1 + ns : 2 * ns, 1 : ns);
       JH = Q' * [sol_update(i).H_til; sol_local(i + 1).H];
       sol_update(i).J = JH(1 : ns, 1 : np);
       sol_update(i + 1).H_til = JH(1 + ns : 2 * ns, 1 : np);
       db = Q' * [sol_update(i).b_til; sol_local(i + 1).b];
       sol_update(i).d = db(1 : ns);
       sol_update(i + 1).b_til = db(1 + ns : 2 * ns);
   end
   sol_update(size_sol).A_til = sol_update(size_sol).G_til;
end

sol_small(num_workers + 1) = struct('A', [], 'C', [], 'H', [], 'b', [], 'B', []);
for i = 1 : num_workers
    sol_worker = cell2mat(sol_update(i));
    sol_small(i).A = sol_worker(end).A_til;
    sol_small(i).C = sol_worker(end).C_til;
    sol_small(i).H = sol_worker(end).H_til;
    sol_small(i).b = sol_worker(end).b_til;
end
sol_small(1).B = sol(1).B;
sol_small(end).B = sol(end).B;
sol_small(end).H = sol(end).H;
sol_small(end).b = sol(end).b;
sol_small = sequential_qr_spmd(sol_small);
sol_small = partition_backward_substitution(sol_small);

spmd
    sol_update(1).delta_s = sol_small(labindex).delta_s;
    sol_update(size_sol + 1).delta_s = sol_small(labindex + 1).delta_s;
    delta_p = sol_small(end).delta_p;
    for i = size_sol : -1 : 2
        sol_update(i).delta_s = sol_update(i - 1).R \ (sol_update(i - 1).d - sol_update(i - 1).G * sol_update(1).delta_s - sol_update(i - 1).E * sol_update(i + 1).delta_s - sol_update(i - 1).J * delta_p);
    end
end
delta_s = [];
for i = 1 : num_workers
   sol_worker = cell2mat(sol_update(i));
   for j = 1 : size(sol_worker, 2) - 1
      delta_s =  [delta_s; sol_worker(j).delta_s];
   end
end
delta_s = [delta_s; sol_small(end).delta_s; sol_small(end).delta_p];
end

function sol = sequential_qr_spmd(sol)

% QR decomposition
global ny nz np

M = size(sol, 2);

ns = ny + nz;

sol(1).C_tilde = sol(1).C;
sol(1).G_tilde = sol(1).A;
sol(1).H_tilde = sol(1).H;
sol(1).b_tilde = sol(1).b;
for i = 1 : M - 2
    [Q, R] = qr([sol(i).C_tilde; sol(i + 1).A]);
    sol(i).R = R(1 : ns, :);
    EC = Q' * [zeros(ns); sol(i + 1).C];
    sol(i).E = EC(1 : ns, 1 : ns);
    sol(i + 1).C_tilde = EC(1 + ns : 2 * ns, 1 : ns);
    GG = Q' * [sol(i).G_tilde; zeros(ns)];
    sol(i).G = GG(1 : ns, 1 : ns);
    sol(i + 1).G_tilde = GG(1 + ns : 2 * ns, 1 : ns);
    JH = Q' * [sol(i).H_tilde; sol(i + 1).H];
    sol(i).J = JH(1 : ns, 1 : np);
    sol(i + 1).H_tilde = JH(1 + ns : 2 * ns, 1 : np);
    db = Q' * [sol(i).b_tilde; sol(i + 1).b];
    sol(i).d = db(1 : ns);
    sol(i + 1).b_tilde = db(1 + ns : 2 * ns);
end
[Q, R] = qr([sol(M - 1).C_tilde sol(M - 1).G_tilde sol(M - 1).H_tilde; sol(M).B sol(1).B sol(M).H]);
sol(M - 1).R = R(1 : ns, 1 : ns);
sol(M - 1).G = R(1 : ns, 1 + ns : 2 * ns);
sol(M - 1).J = R(1 : ns, 1 + 2 * ns : 2 * ns + np);
sol(M).R = R(1 + ns : 2 * ns, 1 + ns : 2 * ns);
sol(M).J = R(1 + ns : 2 * ns, 1 + 2 * ns : 2 * ns + np);
sol(M).Rp = R(1 + 2 * ns : 2 * ns + np, 1 + 2 * ns : 2 * ns + np);

d = Q' * [sol(M - 1).b_tilde; sol(M).b];
sol(M -1).d = d(1 : ns);
sol(M).d = d(1 + ns : 2 * ns);
sol(M).dp = d(1 + 2 * ns : 2 * ns + np);
end

function sol = partition_backward_substitution(sol)

global ny nz 

M = size(sol, 2);

ns = ny + nz;

% delta_p = zeros(np, 1);
delta_p = sol(M).Rp \ sol(M).dp;
sol(M).delta_p = delta_p;
% delta_s1 = zeros(ns, 1);
delta_s1 = sol(M).R \ (sol(M).d - sol(M).J * delta_p);
sol(1).delta_s = delta_s1;
% delta_sN = zeros(ns, 1);
delta_sN = sol(M - 1).R \ (sol(M - 1).d - sol(M-1).G * delta_s1 - sol(M - 1).J * delta_p);
sol(M).delta_s = delta_sN;
for j = M - 1 : -1 : 2
%     delta_sj = zeros(ns, 1);
    delta_sj = sol(j - 1).R \ (sol(j - 1).d - sol(j - 1).E * sol(j + 1).delta_s - sol(j - 1).G * sol(1).delta_s - sol(j - 1).J * sol(M).delta_p);
    sol(j).delta_s = delta_sj;
end
end

%% Mesh refinement of the problem
function [s0, tspan] = mesh_refinement(sol, s0, tspan)
global ny nz M

%% collect all the local truncation error
lte_all = zeros(M - 1, 1);
for i = 1 : M - 1
    lte_all(i) = sol(i).e;
end
%% Deleting Nodes
i = 1;
k_D = 0; % Record the number of the deleted nodes

while i < M - 4
    lte_i = lte_all(i);
    if lte_i <= 1e-3
        lte_iplus2 = lte_all(i + 1);
        lte_iplus3 = lte_all(i + 2);
        lte_iplus4 = lte_all(i + 3);
        lte_iplus5 = lte_all(i + 4);
        if ((lte_iplus2 <= 1e-3) && (lte_iplus3 <= 1e-3) && (lte_iplus4 <= 1e-3) && (lte_iplus5 <= 1e-3))
            tspan(i + 1) = [];
            lte_all(i + 1) = [];
            s0((1 + (i) * (nz + ny) : (i + 1) * (nz + ny))) = [];
            tspan(i + 2) = [];
            lte_all(i+2) = [];
            s0(1 + (i + 1) * (nz + ny) : (i + 2) * (nz + ny)) = [];
            M = M - 2;
            k_D = k_D + 2;
            i = i + 2;
        end
    end
    i = i + 1;
end

%% Adding Nodes

i = 1;
k_A = 0; % Record the number of the added nodes

while i <= M - 1
    lte_i = lte_all(i);
    if lte_i > 1
        if lte_i > 100
            delta_t = (tspan(i + 1) - tspan(i)) / 4;
            t_i = tspan(i);
            t_iplus1 = t_i + delta_t;
            t_iplus2 = t_i + 2 * delta_t;
            t_iplus3 = t_i + 3 * delta_t;
            tspan = [tspan(1 : i), t_iplus1, t_iplus2, t_iplus3, tspan(i+1 : end)];
            if i == M - 1 
                delta_lte = lte_all(i) / 4;
            else
                delta_lte = (lte_all(i + 1) - lte_all(i)) / 4;
            end
            lte_iplus1 = lte_i + delta_lte;
            lte_iplus2 = lte_i + 2 * delta_lte;
            lte_iplus3 = lte_i + 3 * delta_lte;
            lte_all = [lte_all(1 : i); lte_iplus1; lte_iplus2; lte_iplus3; lte_all(i+1:end)];
            delta_s0 = (s0(1 + (i) * (nz + ny) : (i + 1) * (nz + ny)) - s0(1 + (i - 1)*(nz+ny) : (i)*(nz+ny))) / 4;
            s0_i = s0(1 + (i - 1)*(nz+ny) : (i)*(nz+ny));
            s0_iplus1 = s0_i + delta_s0;
            s0_iplus2 = s0_i + 2 * delta_s0;
            s0_iplus3 = s0_i + 3 * delta_s0;
            s0 = [s0(1 : (i)*(nz+ny));  s0_iplus1; s0_iplus2; s0_iplus3; s0(1 + i * (nz + ny):end)];
            M = M + 3;
            k_A = k_A + 3;
            i = i + 3;
        else
            delta_t = (tspan(i + 1) - tspan(i)) / 2;
            t_i = tspan(i);
            t_iplus1 = t_i + delta_t;
            tspan = [tspan(1 : i), t_iplus1, tspan(i+1 : end)];
            if i == M - 1 
                delta_lte = lte_all(i) / 2;
            else
                delta_lte = (lte_all(i + 1) - lte_all(i)) / 2;
            end
            lte_iplus1 = lte_i + delta_lte;
            lte_all = [lte_all(1 : i); lte_iplus1; lte_all(i+1:end)];
            delta_s0 = (s0(1 + (i) * (nz + ny) : (i + 1) * (nz + ny)) - s0(1 + (i - 1)*(nz+ny) : (i)*(nz+ny))) / 2;
            s0_i = s0(1 + (i - 1)*(nz+ny) : (i)*(nz+ny));
            s0_iplus1 = s0_i + delta_s0;
            s0 = [s0(1 : (i)*(nz+ny));  s0_iplus1; s0(1 + i * (nz + ny):end)];
            M = M + 1;
            k_A = k_A + 1;
            i = i + 1;
        end
    end
    i = i + 1;
end
end

%% Convert the solution from matrix form to vector form
function s = matrixtovec(y, z)

global ny nz M

s = ones(M * (ny + nz), 1);
for i = 1 : ny
    s(i : ny + nz : end) = y(:, i);
end
for i = 1 : nz
    s(ny + i : ny + nz : end) = z(:, i);
end
end

%% Convert the solution from vector form to matrix form
function [y, z] = vectomatrix(s0)

global ny nz M

y = zeros(M, ny);
z = zeros(M, nz);
for i = 1 : ny
    y(:, i) = s0(i : ny + nz : end);
end
for i = 1 : nz
    z(:, i) = s0(ny + i : ny + nz : end);
end
end

%% Row method for DAE integration
function [x_next, lte] = row_step(ODE_h, DAE_g, J_hg, tspan, x0, p, alpha0, ny, nz, np, coefficients, tol)
%% Single step integration of the DAE
y0 = x0(1 : ny);
z0 = x0(1 + ny : end);
delta = tspan(end) - tspan(1);
M = zeros(ny + nz);
M(1 : ny, 1 : ny) = eye(ny);
W_j = J_hg(tspan, x0, p, alpha0);
[L,U] = lu(M - coefficients.gamma(1, 1) * delta * W_j);

v_j = zeros(ny + nz , 1);
v_j(ny + 1 : end) = DAE_g(y0, z0, p, alpha0);
x_stage = zeros(ny + nz, coefficients.stages);
k = zeros(ny + nz, coefficients.stages);
f = zeros(ny + nz, coefficients.stages);
x_stage(:, 1) = x0;
f(:, 1) = [ODE_h(y0, z0, p, alpha0); DAE_g(y0, z0, p, alpha0)];
k(:, 1) = U \ (L \ (f(:, 1) - v_j));
for i = 2 : coefficients.stages
    sum1 = zeros(ny + nz, 1);
    sum2 = zeros(ny + nz, 1);
    for l = 1 : i - 1
        sum1 = sum1 + coefficients.alpha(i, l) * k(:, l);
        sum2 = sum2 + coefficients.gamma(i, l) * k(:, l);
    end
    x_stage(:, i) = x0 + delta * sum1;
    y_i = x_stage(1 : ny, i);
    z_i = x_stage(ny + 1 : end, i);
    f(:, i) = [ODE_h(y_i, z_i, p, alpha0); DAE_g(y_i, z_i, p, alpha0)];
    k(:, i) = U \ (L \ (f(:, i) + delta * W_j * sum2 - v_j));
end
sum3 = zeros(ny + nz, 1);
for i =1 : coefficients.stages
    sum3 = sum3 + coefficients.b(i) * k(:, i);
end
x_next = x0 + delta * sum3;

%% Local truncation error of the integration result
sigma = 0;
sum4 = zeros(ny + nz, 1);
for i = 1 : coefficients.stages
    sum4 = sum4 + coefficients.e(i) * k(:, i);
end
error_sum = delta * sum4;
for i = 1 : ny + nz
    fac = abs(error_sum(i)) / (tol * (1.0 + abs(x_next(i))));
    sigma = sigma + fac*fac;
end
lte = sqrt(sigma / (ny + nz));
% lte_all = norm(error_sum, Inf) / (tol * (1 + x_next));
% lte = lte_all / (ny + nz);
end

%% Row method for DAE integration with sensitivity
function [x_next, X_next] = row_sensitivity_step(ODE_h, DAE_g, D_hg, J_hg,tspan, x0, p, alpha0, ny, nz, np, coefficients)
% Single step integration of the DAE
y0 = x0(1 : ny);
z0 = x0(1 + ny : end);
delta = tspan(end) - tspan(1);
M = zeros(ny + nz);
M(1 : ny, 1 : ny) = eye(ny);
W_j = J_hg(tspan, x0, p, alpha0);
[L,U] = lu(M - coefficients.gamma(1, 1) * delta * W_j);

v_j = zeros(ny + nz , 1);
v_j(ny + 1 : end) = DAE_g(y0, z0, p, alpha0);
x_stage = zeros(ny + nz, coefficients.stages);
k = zeros(ny + nz, coefficients.stages);
f = zeros(ny + nz, coefficients.stages);
x_stage(:, 1) = x0;
f(:, 1) = [ODE_h(y0, z0, p, alpha0); DAE_g(y0, z0, p, alpha0)];
k(:, 1) = U \ (L \ (f(:, 1) - v_j));
for i = 2 : coefficients.stages
    sum1 = zeros(ny + nz, 1);
    sum2 = zeros(ny + nz, 1);
    for l = 1 : i - 1
        sum1 = sum1 + coefficients.alpha(i, l) * k(:, l);
        sum2 = sum2 + coefficients.gamma(i, l) * k(:, l);
    end
    x_stage(:, i) = x0 + delta * sum1;
    y_i = x_stage(1 : ny, i);
    z_i = x_stage(ny + 1 : end, i);
    f(:, i) = [ODE_h(y_i, z_i, p, alpha0); DAE_g(y_i, z_i, p, alpha0)];
    k(:, i) = U \ (L \ (f(:, i) + delta * W_j * sum2 - v_j));
end
sum3 = zeros(ny + nz, 1);
for i =1 : coefficients.stages
    sum3 = sum3 + coefficients.b(i) * k(:, i);
end
x_next = x0 + delta * sum3;

% Single step integration of the sensitivity
[~,~,h_p,g_y,g_z,g_p] = D_hg(x0, p, alpha0);
V_j = zeros(ny + nz, ny + nz + np);
V_j(1 + ny : end, :) = [g_y g_z g_p];
X0 = zeros(ny + nz, ny + nz + np);
X0(1  :ny, 1 : ny) = eye(ny);
X0(1 + ny : end, 1 + ny : ny + nz) = eye(nz);
X_stage = zeros(ny + nz, coefficients.stages * (ny + nz + np));
K = zeros(ny + nz, coefficients.stages * (ny + nz + np));
X_stage(:, 1 : (ny + nz + np)) = X0;
F_p = zeros(ny + nz, ny + nz + np);
F_p(:, 1 + ny + nz : end) = [h_p; g_p];
W = jacobian_DAE(tspan, x_stage(:, 1), p, alpha0);
K(:, 1 : (ny + nz + np)) = U \ (L \ (W * X_stage(:, 1 : (ny + nz + np)) + F_p - V_j));
for i = 2 : coefficients.stages
    Sum1 = zeros(ny + nz, ny + nz + np);
    Sum2 = zeros(ny + nz, ny + nz + np);
    for l = 1:i - 1
        Sum1 = Sum1 + coefficients.alpha(i, l) * K(:, 1 + (l - 1) * (ny + nz + np): l * (ny + nz + np));
        Sum2 = Sum2 + coefficients.gamma(i, l) * K(:, 1 + (l - 1) * (ny + nz + np): l * (ny + nz + np));
    end
    X_stage(:, 1 + (i - 1) * (ny + nz + np): i * (ny + nz + np)) = X0 + delta * Sum1;
    [~,~,h_p,~,~,g_p] = D_hg(x_stage(: , i), p, alpha0);
    F_p = zeros(ny + nz, ny + nz + np);
    F_p(:, 1 + ny + nz : end) = [h_p; g_p];
    W = J_hg(tspan, x_stage(:, i), p, alpha0);
    K(:, 1 + (i - 1) * (ny + nz + np): i * (ny + nz + np)) = U \ (L \ (W *  X_stage(:, 1 + (i - 1) * (ny + nz + np): i * (ny + nz + np)) + F_p - V_j + delta * W_j * Sum2));
end
Sum3 = zeros(ny + nz, ny + nz + np);
for i =1 : coefficients.stages
    Sum3 = Sum3 + coefficients.b(i) * K(:, 1 + (i - 1) * (ny + nz + np): i * (ny + nz + np));
end
X_next = X0 + delta * Sum3;
end

%% Sequential QR solver
function sol = sequential_qr(sol)

% QR decomposition
global ny nz np M

ns = ny + nz;

sol(1).C_tilde = sol(1).C;
sol(1).G_tilde = sol(1).A;
sol(1).H_tilde = sol(1).H;
sol(1).b_tilde = sol(1).b;
for i = 1 : M - 2
    [Q, R] = qr([sol(i).C_tilde; sol(i + 1).A]);
    sol(i).R = R(1 : ns, :);
    EC = Q' * [zeros(ns); sol(i + 1).C];
    sol(i).E = EC(1 : ns, 1 : ns);
    sol(i + 1).C_tilde = EC(1 + ns : 2 * ns, 1 : ns);
    GG = Q' * [sol(i).G_tilde; zeros(ns)];
    sol(i).G = GG(1 : ns, 1 : ns);
    sol(i + 1).G_tilde = GG(1 + ns : 2 * ns, 1 : ns);
    JH = Q' * [sol(i).H_tilde; sol(i + 1).H];
    sol(i).J = JH(1 : ns, 1 : np);
    sol(i + 1).H_tilde = JH(1 + ns : 2 * ns, 1 : np);
    db = Q' * [sol(i).b_tilde; sol(i + 1).b];
    sol(i).d = db(1 : ns);
    sol(i + 1).b_tilde = db(1 + ns : 2 * ns);
end
[Q, R] = qr([sol(M - 1).C_tilde sol(M - 1).G_tilde sol(M - 1).H_tilde; sol(M).B sol(1).B sol(M).H]);
sol(M - 1).R = R(1 : ns, 1 : ns);
sol(M - 1).G = R(1 : ns, 1 + ns : 2 * ns);
sol(M - 1).J = R(1 : ns, 1 + 2 * ns : 2 * ns + np);
sol(M).R = R(1 + ns : 2 * ns, 1 + ns : 2 * ns);
sol(M).J = R(1 + ns : 2 * ns, 1 + 2 * ns : 2 * ns + np);
sol(M).Rp = R(1 + 2 * ns : 2 * ns + np, 1 + 2 * ns : 2 * ns + np);

d = Q' * [sol(M - 1).b_tilde; sol(M).b];
sol(M -1).d = d(1 : ns);
sol(M).d = d(1 + ns : 2 * ns);
sol(M).dp = d(1 + 2 * ns : 2 * ns + np);
end

function delta_s = backwardsubstitution(sol)

global ny nz np M 

ns = ny + nz;

% delta_p = zeros(np, 1);
delta_p = sol(M).Rp \ sol(M).dp;
sol(M).delta_p = delta_p;
% delta_s1 = zeros(ns, 1);
delta_s1 = sol(M).R \ (sol(M).d - sol(M).J * delta_p);
sol(1).delta_s = delta_s1;
% delta_sN = zeros(ns, 1);
delta_sN = sol(M - 1).R \ (sol(M - 1).d - sol(M-1).G * delta_s1 - sol(M - 1).J * delta_p);
sol(M).delta_s = delta_sN;
for j = M - 1 : -1 : 2
%     delta_sj = zeros(ns, 1);
    delta_sj = sol(j - 1).R \ (sol(j - 1).d - sol(j - 1).E * sol(j + 1).delta_s - sol(j - 1).G * sol(1).delta_s - sol(j - 1).J * sol(M).delta_p);
    sol(j).delta_s = delta_sj;
end
delta_s = zeros(M * (ny + nz) + np, 1);
for i = 1 : M
    delta_s(1 + (i - 1) * ns : i * ns) = sol(i).delta_s;
end
delta_s(1 + M * (ny + nz) : M * (ny + nz) + np) = sol(M).delta_p;
end