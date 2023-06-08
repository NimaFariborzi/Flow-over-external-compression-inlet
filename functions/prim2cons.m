function U = prim2cons(rho,u,v,T,cv)
%PRIM2CONS Converts primitive variables to state space vector of conservative
%variables
    nx = size(rho,1);
    ny = size(rho,2);
    U = zeros(4,nx,ny);
    U(1,:,:) = rho;
    U(2,:,:) = rho.*u;
    U(3,:,:) = rho.*v;
    U(4,:,:) = rho.*(cv*T + 0.5*(u.^2 + v.^2));
end