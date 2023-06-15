clc; clear; close all

% ======================== %
% GRID LAYOUT              %
% ======================== %
%{
     _______________________
    |       |       |       |
    | Block | Block | Block |
    |   4   |   5   |   6   |
    |_______|_______|_______|
    |       |       |       |
    | Block | Block | Block |
    |   1   |   2   |   3   |
    |_______|_______|_______|

%}
% ======================== %
% GRID PARAMETERS          %
% ======================== %

% Angles
delta1 = 14.9766; % Ramp angle 1
delta2 = 18.8126; % Ramp angle 2

% Distances
x1_dist = 1;             % Length of block 1
x2_x1   = 1.1631;        % Ratio of block 1 length to block 2
x2_dist = x1_dist*x2_x1; % Length of block 2
x3_dist = 1;             % Length of block 3

y1_x1   = 1.2395;        % Ratio of block 1 height to length
y1_dist = x1_dist*y1_x1; % Height of blocks 1,2,3
y2_dist = 1.5;           % Height of blocks 4,5

% Computational grid scale
gs = 6.3e-6;

x_int = 1.9670*x1_dist*gs; % x-location of shock intersection  

% Grid spacing
% d_xi = 0.01 * gs; % Grid spacing in xi
d_xi = 0.0125 * gs;
d_et = 0.0125 * gs; % Grid spacing in eta

% Grid stretching
alpha1 = 0.85;
alpha2 = 0.9;

% ======================== %
% BLOCK CONSTRUCTION       %
% ======================== %

% Block one: computational grid
xi1 = 0:d_xi:1*gs;
et1 = 0:d_et:1.25*gs;
[XI1, ET1] = ndgrid(xi1,et1);

% Block one: physical grid
X1 = XI1*x1_dist;
Y1 = zeros(size(ET1));

y_max = y1_dist*gs;
for i = 1:length(xi1)
    x = X1(i,1);
    y_min = x*tand(delta1);
    Y1(i,:) = TwoWayBiasY(Y1(i,:),y_min,y_max,alpha1);
end

% Block two: computational grid
xi2 = 1*gs:d_xi:2.2*gs;
et2 = 0:d_et:1.25*gs; % Brutus?
[XI2, ET2] = ndgrid(xi2,et2);

% Block two: physical grid
X2 = rescale(XI2,x1_dist*gs,(x1_dist + x2_dist)*gs);
Y2 = zeros(size(ET2));

x_min = min(X2,[],"all");
x_max = max(X2,[],"all");
for i = 1:length(et2)
    X2(:,i) = OneWayBiasX(X2(:,i),x_min,x_max);
end

y_max = y1_dist*gs;
for i = 1:length(xi2)
    x = X2(i,1);
    y_min = x1_dist*gs*tand(delta1) + (x - x1_dist*gs)*tand(delta2);
    Y2(i,:) = TwoWayBiasY(Y2(i,:),y_min,y_max,alpha1);
end

% Block three: computational grid
xi3 = 2.2*gs:d_xi:3.2*gs;
et3 = 0:d_et:1.25*gs;
[XI3, ET3] = ndgrid(xi3,et3);

% Block three: physical grid
X3 = rescale(XI3,(x1_dist + x2_dist)*gs,(x1_dist + x2_dist + x3_dist)*gs);
Y3 = zeros(size(ET3));

y_min = Y2(end,1);
y_max = y1_dist*gs;
for i = 1:length(xi3)
    Y3(i,:) = TwoWayBiasY(Y3(i,:),y_min,y_max,alpha1);
end

% Block four: computational grid
xi4 = 0:d_xi:1*gs;
et4 = 1.25*gs:d_et:2.75*gs;
[XI4, ET4] = ndgrid(xi4,et4);

% Block four: physical grid
X4 = XI4*x1_dist;
Y4 = zeros(size(ET4));

y_min = y1_dist*gs;
y_max = (y1_dist + y2_dist)*gs;
for i = 1:length(xi4)
    Y4(i,:) = OneWayBiasY(Y4(i,:),y_min,y_max,alpha2);
end

% Block five: computational grid
xi5 = 1*gs:d_xi:2.2*gs;
et5 = 1.25*gs:d_et:2.75*gs;
[XI5, ET5] = ndgrid(xi5,et5);

% Block five: physical grid
X5 = rescale(XI5,x1_dist*gs,(x1_dist + x2_dist)*gs);
Y5 = zeros(size(ET5));

x_min = min(X5,[],"all");
x_max = max(X5,[],"all");
for i = 1:length(et5)
    X5(:,i) = OneWayBiasX(X5(:,i),x_min,x_max);
end

y_min = y1_dist*gs;
y_max = (y1_dist + y2_dist)*gs;
for i = 1:length(xi5)
    Y5(i,:) = OneWayBiasY(Y5(i,:),y_min,y_max,alpha2);
end

% Block six: computational grid
xi6 = 2.2*gs:d_xi:3.2*gs;
et6 = 1.25*gs:d_et:2.75*gs;
[XI6, ET6] = ndgrid(xi6,et6);

% Block six: physical grid
X6 = rescale(XI6,(x1_dist + x2_dist)*gs,(x1_dist + x2_dist + x3_dist)*gs);
Y6 = zeros(size(ET6));

y_min = y1_dist*gs;
y_max = (y1_dist + y2_dist)*gs;
for i = 1:length(xi6)
    Y6(i,:) = OneWayBiasY(Y6(i,:),y_min,y_max,alpha2);
end

% Concatenate (remove duplicate points!)
XI = [XI1,          XI4(:,2:end);
      XI2(2:end,:), XI5(2:end,2:end);
      XI3(2:end,:), XI6(2:end,2:end)];
ET = [ET1,          ET4(:,2:end);
      ET2(2:end,:), ET5(2:end,2:end);
      ET3(2:end,:), ET6(2:end,2:end)];
X =  [X1,           X4(:,2:end);
      X2(2:end,:),  X5(2:end,2:end);
      X3(2:end,:),  X6(2:end,2:end)];
Y =  [Y1,           Y4(:,2:end);
      Y2(2:end,:),  Y5(2:end,2:end);
      Y3(2:end,:),  Y6(2:end,2:end)];

% Get row and col num of inlet cowl
cowl_cols = [size(X1,2)-1, size(X1,2)];
x = X(:,1);
cowl_rows = find(x>x_int,1):size(X,1);

% Plot
tiledlayout(1,2)
% Computational grid
nexttile
hold on
plot(XI1(:),ET1(:),'.','MarkerSize',15)
plot(XI2(:),ET2(:),'.','MarkerSize',15)
plot(XI3(:),ET3(:),'.','MarkerSize',15)
plot(XI4(:),ET4(:),'.','MarkerSize',15)
plot(XI5(:),ET5(:),'.','MarkerSize',15)
plot(XI6(:),ET6(:),'.','MarkerSize',15)
% plot(XI(:),ET(:),'.','MarkerSize',15)
title('Computational Grid')
xlabel('$\xi$')
ylabel('$\eta$')
axis equal tight
% Physical grid
nexttile
hold on
plot(X1(:),Y1(:),'.','MarkerSize',15)
plot(X2(:),Y2(:),'.','MarkerSize',15)
plot(X3(:),Y3(:),'.','MarkerSize',15)
plot(X4(:),Y4(:),'.','MarkerSize',15)
plot(X5(:),Y5(:),'.','MarkerSize',15)
plot(X6(:),Y6(:),'.','MarkerSize',15)
% Plot cowl
plot(X(cowl_rows,cowl_cols),Y(cowl_rows,cowl_cols),'k.','MarkerSize',15)
% plot(X(:),Y(:),'.','MarkerSize',15)
title('Physical Grid')
xlabel('$x$')
ylabel('$y$')
axis equal tight

exportgraphics(gcf,'Mesh.png','Resolution',300)

nx = size(X,1);
ny = size(X,2);
save('mesh3_fine.mat','XI','ET','X','Y','d_xi','d_et','nx','ny',...
     'cowl_cols','cowl_rows')

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

function x = OneWayBiasX(x,x_min,x_max)
    n = length(x);
    j = (1:n)*2*pi;
    % This is totally ad hoc but works great for this particular mesh...
    x = sin(j/n-1.55) + 1.8*j/n;
    x = rescale(x,x_min,x_max);
end