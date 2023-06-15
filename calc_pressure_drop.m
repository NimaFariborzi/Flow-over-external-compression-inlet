clc; clear; close all

% Load mesh and results
load('mesh3.mat')
load('results.mat')

% Define parameters
M_inf = 3;
T_inf = 288.15; % K, STP
P_inf = 101.3e3; % Pa, STP
cp = 1005; % J/(kg K)
cv = 718; % J/(kg K)
R = cp - cv; % J/(kg K)
gamma = cp/cv; % Unitless, heat capacity ratio
a_inf = sqrt(gamma*R*T_inf);
u_inf = M_inf*a_inf;
rho_inf = P_inf/(R*T_inf);

% Get density
rho = P./(R*T);

% Get total pressure of incoming air
P0_in = P_inf*(1 + (gamma-1)/2*M_inf^2) ^ (gamma/(gamma-1))

% Calculat stagnation pressure everywhere
for i = 1:size(X,1)
    for j = 1:size(X,2)
        v_tot2 = u(i,j)^2 + v(i,j)^2;
        M2 = v_tot2/(gamma*R*T(i,j));
        P0(i,j) = P(i,j)*(1 + (gamma-1)/2*M2) ^ (gamma/(gamma-1));
    end
end

figure()
pcolor(X,Y,P0)
colorbar
axis equal tight
xlabel('$x$')
ylabel('$y$')
title('Stagnation Pressure $P_0$')
% exportgraphics(gcf,'Stag_Press.png','Resolution',300)

% Get total pressure of air along outlet inside intake
n = cowl_cols(1);
i = size(X,1);
P0_out  = zeros(1,n);
rho_out = zeros(1,n);
y_out   = zeros(1,n);
for j = 1:n
    P0_out(j) = P0(i,j);
    rho_out(j) = rho(i,j);
    y_out(j) = Y(i,j);
end

% Calculate mass-averaged stagnation pressure at outlet
P0_out = trapz(y_out,rho_out.*P0_out)/trapz(y_out,rho_out)

% Get pressure recovery
P_ratio = P0_out/P0_in