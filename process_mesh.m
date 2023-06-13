clc; clear; close all

% Load mesh and functions
load('mesh3.mat')
addpath('./functions/')

% Calculate inverse metrics
x_xi = ddxi_central(X,d_xi);
y_xi = ddxi_central(Y,d_xi);
x_et = ddet_central(X,d_et);
y_et = ddet_central(Y,d_et);

% Use one-sided differences in eta at cowl top/bottom
% NB: Only actually need to correct y_et (x_et is zero either way)
j = cowl_cols(1); % Underside of cowl
for i = cowl_rows
    % Backward difference
    y_et(i,j) = (3*Y(i,j) - 4*Y(i,j-1) + Y(i,j-2))/(2*d_et);
end
j = cowl_cols(2); % Topside of cowl
for i = cowl_rows
    % Forward difference
    y_et(i,j) = (-3*Y(i,j) + 4*Y(i,j+1) - Y(i,j+2))/(2*d_et);
end

% Calcualte Jacobian
J = x_xi.*y_et - x_et.*y_xi;

% Calculate metrics
xi_x =  y_et./J;
et_x = -y_xi./J;
xi_y = -x_et./J;
et_y =  x_xi./J;

% Export results
save('mesh3_metrics.mat','x_xi','y_xi','x_et','y_et','xi_x','et_x','xi_y','et_y','J')