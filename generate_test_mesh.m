clc; clear; close all

% ======================== %
% GRID PARAMETERS          %
% ======================== %

% Distances
x1_dist = 1e-5; % Length of domain
y1_dist = 8e-6; % Height of domain

% Grid spacing
nx = 75;
ny = 80;

% Grid stretching
alpha_y = 0.7;

% ======================== %
% BLOCK CONSTRUCTION       %
% ======================== %

% computational grid
xi = linspace(0,1e-5,nx);
et = linspace(0,8e-6,ny);
[XI, ET] = ndgrid(xi,et);

d_xi = xi(2) - xi(1);
d_et = et(2) - et(1);

% physical grid
X = XI;
Y = zeros(size(ET));

for i = 1:nx
    Y(i,:) = OneWayBiasY(Y(i,:),0,y1_dist,alpha_y);
end

% Plot
tiledlayout(1,2)
% Computational grid
nexttile
plot(XI(:),ET(:),'.','MarkerSize',15)
title('Computational Grid')
xlabel('$\xi$')
ylabel('$\eta$')
axis equal tight
% Physical grid
nexttile
plot(X(:),Y(:),'.','MarkerSize',15)
title('Physical Grid')
xlabel('$x$')
ylabel('$y$')
axis equal tight

save('test_mesh.mat','XI','ET','X','Y','d_xi','d_et','nx','ny')

function y = TwoWayBiasY(y,y_min,y_max,alpha)
    n = length(y);
    j = 1:n;
    e = -1+2*j/n;
    y = 1/alpha * tanh(e * atanh(alpha)) + 1;
    y = rescale(y,y_min,y_max);
end

function y = OneWayBiasY(y,y_min,y_max,alpha)
    n = length(y);
    j = 1:n;
    e = -1+j/n;
    y = 1/alpha * tanh(e * atanh(alpha)) + 1;
    y = rescale(y,y_min,y_max);
end