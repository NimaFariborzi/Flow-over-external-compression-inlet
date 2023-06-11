% Function to apply BC's to primitive variables
function [u, v, P, T, U] = apply_BCs_flat(u, v, P, T, R, cv, u_inf, P_inf, T_inf, ...
                                     AdiabaticWallFlag)
    % Inlet
    u(1,:) = u_inf;
    v(1,:) = 0; % Accounts for leading edge as well
    P(1,:) = P_inf; % Accounts for leading edge as well
    T(1,:) = T_inf; % Accounts for leading edge as well
    % Far-field
    u(:,end) = u_inf;
    v(:,end) = 0;
    P(:,end) = P_inf;
    T(:,end) = T_inf;
    % Outlet
    u(end,2:end-1) = 2*u(end-1,2:end-1) - u(end-2,2:end-1);
    v(end,2:end-1) = 2*v(end-1,2:end-1) - v(end-2,2:end-1);
    P(end,2:end-1) = 2*P(end-1,2:end-1) - P(end-2,2:end-1);
    T(end,2:end-1) = 2*T(end-1,2:end-1) - T(end-2,2:end-1);
    % Wall
    u(:,1) = 0; % Accounts for leading edge as well
    v(2:end,1) = 0;
    P(2:end,1) = 2*P(2:end,2) - P(2:end,3);
    if ~AdiabaticWallFlag
        T(2:end,1) = T_inf;
    else
        T(2:end,1) = T(2:end,2);
    end
    % Apply changes to the primitive variables to the conservatives
    rho = P./(R*T);
    U = prim2cons(rho,u,v,T,cv);
end