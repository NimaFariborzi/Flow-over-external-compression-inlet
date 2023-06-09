%% Midterm Part 1: Supersonic Flow over Flat Plate
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

%% Functions

% Function to calculate the values of E and F to be used during predictor step
% calculations
function [E, F] = calc_EF_pred(U,u,v,P,T,cp,Pr,dx,dy)
    % Allocate E and F
    E = zeros(size(U));
    F = zeros(size(U));
    % Calculate viscosity and thermal conductivity
    mu = sutherland(T);
    k = mu*cp/Pr;
    % Calculate normal stresses
    tau_xx = mu.*(4/3*ddxi_bwd(u,dx) - 2/3*ddet_central(v,dy));
    tau_yy = mu.*(4/3*ddet_bwd(v,dy) - 2/3*ddxi_central(u,dx));
    % Calculate shear stresses (different for E and F due to finite differences)
    tau_xy_E = mu.*(ddet_central(u,dy) + ddxi_bwd(v,dx));
    tau_xy_F = mu.*(ddxi_central(v,dx) + ddet_bwd(u,dy));
    % Calculate heat fluxes
    q_x = -k.*ddxi_bwd(T,dx);
    q_y = -k.*ddet_bwd(T,dy);
    % First slice
    E(1,:,:) = U(2,:,:);
    F(1,:,:) = U(3,:,:);
    % Second slice
    E(2,:,:) = squeeze(U(2,:,:)).*u + P - tau_xx;
    F(2,:,:) = squeeze(U(2,:,:)).*v     - tau_xy_F;
    % Third slice
    E(3,:,:) = squeeze(U(3,:,:)).*u     - tau_xy_E;
    F(3,:,:) = squeeze(U(3,:,:)).*v + P - tau_yy;
    % Fourth slice
    E(4,:,:) = (squeeze(U(4,:,:)) + P).*u - u.*tau_xx - v.*tau_xy_E + q_x;
    F(4,:,:) = (squeeze(U(4,:,:)) + P).*v - v.*tau_yy - u.*tau_xy_F + q_y;
end

% Function to calculate the values of E and F to be used during corrector step
% calculations
function [E, F] = calc_EF_corr(U,u,v,P,T,cp,Pr,dx,dy)
    % Allocate E and F
    E = zeros(size(U));
    F = zeros(size(U));
    % Calculate viscosity and thermal conductivity
    mu = sutherland(T);
    k = mu*cp/Pr;
    % Calculate normal stresses
    tau_xx = mu.*(4/3*ddxi_fwd(u,dx) - 2/3*ddet_central(v,dy));
    tau_yy = mu.*(4/3*ddet_fwd(v,dy) - 2/3*ddxi_central(u,dx));
    % Calculate shear stresses (different for E and F due to finite differences)
    tau_xy_E = mu.*(ddet_central(u,dy) + ddxi_fwd(v,dx));
    tau_xy_F = mu.*(ddxi_central(v,dx) + ddet_fwd(u,dy));
    % Calculate heat fluxes
    q_x = -k.*ddxi_fwd(T,dx);
    q_y = -k.*ddet_fwd(T,dy);
    % First slice
    E(1,:,:) = U(2,:,:);
    F(1,:,:) = U(3,:,:);
    % Second slice
    E(2,:,:) = squeeze(U(2,:,:)).*u + P - tau_xx;
    F(2,:,:) = squeeze(U(2,:,:)).*v     - tau_xy_F;
    % Third slice
    E(3,:,:) = squeeze(U(3,:,:)).*u     - tau_xy_E;
    F(3,:,:) = squeeze(U(3,:,:)).*v + P - tau_yy;
    % Fourth slice
    E(4,:,:) = (squeeze(U(4,:,:)) + P).*u - u.*tau_xx - v.*tau_xy_E + q_x;
    F(4,:,:) = (squeeze(U(4,:,:)) + P).*v - v.*tau_yy - u.*tau_xy_F + q_y;
end

% Function to apply BC's to primitive variables
function [u, v, P, T, U] = apply_BCs(u, v, P, T, R, cv, u_inf, P_inf, T_inf, ...
                                     AdiabaticWallFlag)
    % Inlet
    u(1,:) = u_inf;
    v(1,:) = 0; % Accounts for leading edge as well
    P(1,:) = P_inf; % Accounts for leading edge as well
    T(1,:) = T_inf; % Accounts for leading edge as well
    % Far-field
    u(:,end) = u_inf;
    v(:,end) = 0;
    P(:,end) = P_inf;
    T(:,end) = T_inf;
    % Outlet
    u(end,2:end-1) = 2*u(end-1,2:end-1) - u(end-2,2:end-1);
    v(end,2:end-1) = 2*v(end-1,2:end-1) - v(end-2,2:end-1);
    P(end,2:end-1) = 2*P(end-1,2:end-1) - P(end-2,2:end-1);
    T(end,2:end-1) = 2*T(end-1,2:end-1) - T(end-2,2:end-1);
    % Wall
    u(:,1) = 0; % Accounts for leading edge as well
    v(2:end,1) = 0;
    P(2:end,1) = 2*P(2:end,2) - P(2:end,3);
    if ~AdiabaticWallFlag
        T(2:end,1) = T_inf;
    else
        T(2:end,1) = T(2:end,2);
    end
    % Apply changes to the primitive variables to the conservatives
    rho = P./(R*T);
    U = prim2cons(rho,u,v,T,cv);
end