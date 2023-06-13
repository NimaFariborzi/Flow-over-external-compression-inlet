clc; clear; close all

% Load mesh and functions
load('mesh3.mat')
addpath('./functions/')

% Calculate inverse metrics
x_xi = ddxi_central(X,d_xi);
y_xi = ddxi_central(Y,d_xi);
x_et = ddet_central(X,d_et);
y_et = ddet_central(Y,d_et,cowl_rows,cowl_cols);
% NB: Only actually need to correct y_et with one-sided differences at the cowl
% since x_et is zero either way, and xi won't have any difference

% Calcualte Jacobian
J = x_xi.*y_et - x_et.*y_xi;

% Calculate metrics
xi_x =  y_et./J;
et_x = -y_xi./J;
xi_y = -x_et./J;
et_y =  x_xi./J;

% Export results
save('mesh3_metrics.mat','x_xi','y_xi','x_et','y_et','xi_x','et_x','xi_y','et_y','J')