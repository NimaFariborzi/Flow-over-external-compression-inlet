clc; clear; close all

% Load results and mesh
load('results.mat')
load('mesh3_fine.mat')

% % x-velocity
% figure()
% pcolor(X,Y,u);
% title('x-velocity')
% cb = colorbar;
% cb.Label.Interpreter = 'latex';
% ylabel(cb,'$u$ $\left[\frac{m}{s} \right]$')
% xlabel('$x$ (m)')
% ylabel('$y$ (m)')
% axis equal tight
% pause(0.2)
% exportgraphics(gcf,'fine_u.png','Resolution',300)
% close
% 
% % y-velocity
% figure()
% pcolor(X,Y,v);
% title('y-velocity')
% cb = colorbar;
% cb.Label.Interpreter = 'latex';
% ylabel(cb,'$v$ $\left[\frac{m}{s} \right]$')
% xlabel('$x$ (m)')
% ylabel('$y$ (m)')
% axis equal tight
% pause(0.2)
% exportgraphics(gcf,'fine_v.png','Resolution',300)
% close

% density
% figure()
% pcolor(X,Y,squeeze(U(1,:,:)));
% title('Density')
% cb = colorbar;
% cb.Label.Interpreter = 'latex';
% ylabel(cb,'$\rho$ $\left[\frac{kg}{m^3} \right]$')
% xlabel('$x$ (m)')
% ylabel('$y$ (m)')
% axis equal tight
% clim([min(squeeze(U(1,:,:)),[],'all'), 10])
% pause(0.2)
% exportgraphics(gcf,'fine_rho.png','Resolution',300)
% close

% pressure
% figure()
% pcolor(X,Y,P);
% title('Pressure')
% cb = colorbar;
% cb.Label.Interpreter = 'latex';
% ylabel(cb,'$P$ $\left[\frac{N}{m^2} \right]$')
% xlabel('$x$ (m)')
% ylabel('$y$ (m)')
% axis equal tight
% clim([min(P,[],'all'), 2e6])
% pause(0.2)
% exportgraphics(gcf,'fine_P.png','Resolution',300)
% close

% temperature
figure()
pcolor(X,Y,T);
title('Temperature')
cb = colorbar;
cb.Label.Interpreter = 'latex';
ylabel(cb,'$T$ $\left[K\right]$')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
axis equal tight
pause(0.2)
exportgraphics(gcf,'fine_T.png','Resolution',300)
close