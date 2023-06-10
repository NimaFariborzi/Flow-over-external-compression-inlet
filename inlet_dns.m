%% Final Project Main File
clc; clear; close all

% Use my functions
addpath('./functions')

% Define parameters
M_inf = 4;
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
dt = 2.35e-11;
num_steps = 1500;

% Create grid
L = 1e-5; % m, length of domain
H = 8e-6; % m, height of domain
nx = 75; % number of grid points in x
ny = 80; % number of grid points in y
x = linspace(0,L,nx); % x vector
y = linspace(0,H,ny); % y vector
dx = x(2) - x(1); % grid spacing in x
dy = y(2) - y(1); % grid spacing in y
[X,Y] = ndgrid(x,y); % create ndgrid
% NB: All ddx/ddy functions are set up to be used with ndgrid arrays, not
% meshgrid

% Allocate solution arrays
u = zeros(nx,ny); % x-velocity
v = zeros(nx,ny); % y-velocity
P = zeros(nx,ny); % pressure
T = zeros(nx,ny); % temperature

% Set initial conditions via primitives
u(:,:) = u_inf;
P(:,:) = P_inf;
T(:,:) = T_inf;

% Apply BC's
[u, v, P, T, U] = apply_BCs(u, v, P, T, R, cv, u_inf, P_inf, T_inf, ...
                            AdiabaticWallFlag);
[~,~,~,~,~,e,~] = cons2prim(U,R,cv); % Get e for plotting

% Visualization parameters
axis_FS = 20;
cb_FS = 25;
title_FS = 30;
plot_frequency = 25;

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
semilogy(1,'o','MarkerEdgeColor','none')
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
    [E, F] = calc_EF_pred(U,u,v,P,T,cp,Pr,dx,dy);
    % Calculate U_pred
    U_pred = U + dt*(-ddxi_fwd_3(E,dx) - ddet_fwd_3(F,dy));
    % Get required primitive variables back from U_pred
    [~, u_pred, v_pred, T_pred, P_pred, ~, ~] = cons2prim(U_pred,R,cv);
    % Apply BC's to all predictor variables
    [u_pred, v_pred, P_pred, T_pred, U_pred] = ...
        apply_BCs(u_pred, v_pred, P_pred, T_pred, R, cv, u_inf, P_inf, ...
                  T_inf, AdiabaticWallFlag);
    % ======================= Calculate corrector step ======================= %
    % Get value of E and F
    [E_pred, F_pred] = ...
        calc_EF_corr(U_pred,u_pred,v_pred,P_pred,T_pred,cp,Pr,dx,dy);
    % Calculate U
    U = 0.5*(U + U_pred) + dt/2*(-ddxi_bwd_3(E_pred,dx) - ddet_bwd_3(F_pred,dy));
    % Get required primitive variables back from U
    [~, u, v, T, P, ~, ~] = cons2prim(U,R,cv);
    % Apply BC's to all variables
    [u, v, P, T, U] = apply_BCs(u, v, P, T, R, cv, u_inf, P_inf, T_inf, ...
                                AdiabaticWallFlag);
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