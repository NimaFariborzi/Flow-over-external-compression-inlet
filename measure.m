clc; clear; close all

% Original ramp for M = 2
% delta1 = 10.405;
% beta1 = 39.758;
% 
% delta2 = 21.539;
% beta2 = 51.898;

% Modified ramp for M = 3
delta1 = 14.9766;
beta1 = 32.2162;

delta2 = 18.8126;
beta2 = 45.1465;

% Modified ramp for M = 4
% delta1 = 16.1345;
% beta1 = 28.2389;
% 
% delta2 = 22.1015;
% beta2 = 41.6599;

alpha1 = beta1 - delta1;
alpha2 = beta2 - delta2;

gamma1 = 180 - beta2 + delta1;
theta1 = 180 - alpha1 - gamma1;

x1 = 1;

a1 = x1/cosd(delta1);
a3 = sind(alpha1)*a1/sind(theta1);

b1 = a3*cosd(alpha2);
x2 = b1*cosd(delta2)

y_int = a1*sind(delta1) + a3*sind(beta2)

% x value of shock intersection
x_int = x1 + a3*cosd(beta2)