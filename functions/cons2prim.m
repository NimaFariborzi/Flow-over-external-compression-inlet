function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
%CONS2PRIM Converts state space vector of conservative variables to primitive
%variables
    rho = squeeze(U(1,:,:));
    u = squeeze(U(2,:,:))./rho;
    v = squeeze(U(3,:,:))./rho;
    Et = squeeze(U(4,:,:));
    e = Et./rho - 0.5*(u.^2 + v.^2);
    T = e/cv;
    p = rho.*(R*T);
end