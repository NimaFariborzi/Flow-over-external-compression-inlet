clc; clear; close all

delta1 = 10.405;
beta1 = 39.758;
alpha1 = beta1 - delta1;

delta2 = 21.539;
beta2 = 51.898;
alpha2 = beta2 - delta2;

gamma1 = 180 - beta2 + delta1;
theta1 = 180 - alpha1 - gamma1;

x1 = 1;

a1 = x1/cosd(delta1);
a3 = sind(alpha1)*a1/sind(theta1);

b1 = a3*cosd(alpha2);
x2 = b1*cosd(delta2)

y = a1*sind(delta1) + a3*sind(beta2)