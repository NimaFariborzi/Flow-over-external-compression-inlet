%% Final Project Main File
clc; clear; close all

% Load functions
addpath('./functions/')

% Load mesh and metrics
load('./mesh3.mat')
load('./mesh3_metrics.mat')

% Define parameters
M_inf = 3;
T_inf = 288.15; % K, STP
P_inf = 101.3e3; % Pa, STP
cp = 1005; % J/(kg K)
cv = 718; % J/(kg K)
R = cp - cv; % J/(kg K)
gamma = cp/cv; % Unitless, heat capacity ratio
Pr = 0.71; % Unitless, Prandtl number
a_inf = sqrt(gamma*R*T_inf);
u_inf = M_inf*a_inf;

% Adiabatic wall flag
AdiabaticWallFlag = false;

% Time controls
t = 0;
dt = 3e-12;
num_steps = 9000;

% Allocate solution arrays
u = zeros(nx,ny); % x-velocity
v = zeros(nx,ny); % y-velocity
P = zeros(nx,ny); % pressure
T = zeros(nx,ny); % temperature

% Allocate intermediate arrays
Eps      = zeros(4,nx,ny);
Phi      = zeros(4,nx,ny);
Upsilon  = zeros(4,nx,ny);
U_pred   = zeros(4,nx,ny);
Eps_pred = zeros(4,nx,ny);
Phi_pred = zeros(4,nx,ny);

% Set initial conditions via primitives
u(:,:) = u_inf;
P(:,:) = P_inf;
T(:,:) = T_inf;

% Apply BC's
[u, v, P, T, U] = apply_BCs_inlet(u, v, P, T, R, cv, u_inf, P_inf, T_inf, X, Y, cowl_rows, cowl_cols);
[~,~,~,~,~,e,~] = cons2prim(U,R,cv); % Get e for plotting

% Optionally, if results from a previous run exist, load them...
load('results.mat')

% Visualization parameters
axis_FS = 20;
cb_FS = 25;
title_FS = 30;
plot_frequency = 100;

% Create solution visualizations
figure(1)
ti_la = tiledlayout(2,3,'TileSpacing','compact');
title(ti_la,'Flow Properties','Interpreter','latex','FontSize',title_FS)
xlabel(ti_la,'$x$ (m)','Interpreter','latex','FontSize',title_FS)
ylabel(ti_la,{'$y$ (m)';' '},'Interpreter','latex','FontSize',title_FS)
time_string = sprintf('$t = $ \\texttt{%0.3e}',t);
time_subtitle = subtitle(ti_la,time_string,'Interpreter','latex', ...
                         'FontSize',title_FS);

% rho
nexttile
rho_plot = pcolor(X,Y,squeeze(U(1,:,:)));
title('Density')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$\rho$ $\left[\frac{kg}{m^3} \right]$')
% clim([1.0738 3.8304]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% u
nexttile
u_plot = pcolor(X,Y,u);
title('x-velocity')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$u$ $\left[\frac{m}{s} \right]$')
% clim([0, u_inf]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% v
nexttile
v_plot = pcolor(X,Y,v);
title('y-velocity')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$v$ $\left[\frac{m}{s} \right]$')
% clim([-0.4479 170.8420]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% e
nexttile
e_plot = pcolor(X,Y,e);
title('Internal Energy')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$e$ $\left[\frac{J}{kg} \right]$')
% clim([1.9606e+05 3.3003e+05]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% P
nexttile
P_plot = pcolor(X,Y,P);
title('Pressure')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$P$ $\left[\frac{N}{m^2} \right]$')
% clim([1e5, 3e5]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% T
nexttile
T_plot = pcolor(X,Y,T);
title('Temperature')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cb_FS;
ylabel(cb,'$T$ $\left[K\right]$')
% clim([273.0592 459.6521]);
ax = gca;
ax.FontSize = axis_FS;
axis equal tight

% Create convergence visualization
f2 = figure(2);
f2.WindowState = 'normal';
% Make invisible point to force axes to semilogy
semilogy(1e-1,'o','MarkerEdgeColor','none')
xlabel('Step')
ylabel('Relative Difference in $\rho$')
title('Convergence')
conv_line = animatedline;

drawnow
U_old = U;
for i = 1:num_steps
    % Increment time
    t = t + dt;
    % ======================= Calculate predictor step ======================= %
    % Get value of E and F
    [E, F] = calc_EF_pred(U,u,v,P,T,cp,Pr,d_xi,d_et,xi_x,et_x,et_y,xi_y,cowl_rows,cowl_cols);
    % Get script E, F, and U
    for j = 1:4
        Eps(j,:,:) =  y_et.*squeeze(E(j,:,:)) - x_et.*squeeze(F(j,:,:));
        Phi(j,:,:) = -y_xi.*squeeze(E(j,:,:)) + x_xi.*squeeze(F(j,:,:));
        Upsilon(j,:,:) = J.*squeeze(U(j,:,:));
    end
    % Calculate Upsilon_pred
    Upsilon_pred = Upsilon + dt*(-ddxi_fwd_3(Eps,d_xi,cowl_rows,cowl_cols) - ddet_fwd_3(Phi,d_et,cowl_rows,cowl_cols));
    % Get U_pred from Upsilon_pred
    for j = 1:4
        U_pred(j,:,:) = squeeze(Upsilon_pred(j,:,:))./J;
    end
    % Get required primitive variables back from U_pred
    [~, u_pred, v_pred, T_pred, P_pred, ~, ~] = cons2prim(U_pred,R,cv);
    % Apply BC's to all predictor variables
    [u_pred, v_pred, P_pred, T_pred, U_pred] = ...
        apply_BCs_inlet(u_pred, v_pred, P_pred, T_pred, R, cv, u_inf, P_inf, ...
                  T_inf, X, Y, cowl_rows, cowl_cols);
    % Update Upsilon_pred to reflect applied BCs
    for j = 1:4
        Upsilon_pred(j,:,:) = J.*squeeze(U_pred(j,:,:));
    end
    % ======================= Calculate corrector step ======================= %
    % Get value of E and F
    [E_pred, F_pred] = ...
        calc_EF_corr(U_pred,u_pred,v_pred,P_pred,T_pred,cp,Pr,d_xi,d_et,xi_x,et_x,et_y,xi_y,cowl_rows,cowl_cols);
    % Get script E_pred and F_pred
    for j = 1:4
        Eps_pred(j,:,:) =  y_et.*squeeze(E_pred(j,:,:)) - x_et.*squeeze(F_pred(j,:,:));
        Phi_pred(j,:,:) = -y_xi.*squeeze(E_pred(j,:,:)) + x_xi.*squeeze(F_pred(j,:,:));
    end
    % Calculate Upsilon
    Upsilon = 0.5*(Upsilon + Upsilon_pred) + dt/2*(-ddxi_bwd_3(Eps_pred,d_xi,cowl_rows,cowl_cols) - ddet_bwd_3(Phi_pred,d_et,cowl_rows,cowl_cols));
    % Get U from Upsilon
    for j = 1:4
        U(j,:,:) = squeeze(Upsilon(j,:,:))./J;
    end
    % Get required primitive variables back from U
    [~, u, v, T, P, ~, ~] = cons2prim(U,R,cv);
    % Apply BC's to all variables
    [u, v, P, T, U] = apply_BCs_inlet(u, v, P, T, R, cv, u_inf, P_inf, T_inf, X, Y, cowl_rows, cowl_cols);
    % ================================= Done! ================================ %
    % Compute convergence via max relative change in density
    max_diff_rel = max(abs((U(1,:,:) - U_old(1,:,:))./U(1,:,:)),[],'all');
    addpoints(conv_line,i,max_diff_rel)
    % Update plots
    if ~mod(i,plot_frequency) || i == num_steps
        [~,~,~,~,~,e,~] = cons2prim(U,R,cv); % Get e for plotting
        rho_plot.CData = squeeze(U(1,:,:));
        u_plot.CData = u;
        v_plot.CData = v;
        e_plot.CData = e;
        P_plot.CData = P;
        T_plot.CData = T;
        % Update subtitle
        time_string = sprintf('$t = $ \\texttt{%0.3e}',t);
        time_subtitle.String = time_string;
        % Display
        drawnow
    end
    U_old = U;
end

save('results.mat','U','u','v','P','T','t')