% Function to apply BC's to primitive variables
function [u, v, P, T, U] = apply_BCs_inlet(u, v, P, T, R, cv, u_inf, P_inf, T_inf, ...
                                     AdiabaticWallFlag)

    % Inlet and far-field are same as they were for flat plate, since they
    % are Dirichlet BC's

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

    % Outlet can remain the same as well since there is no grid-stretching
    % in the x-direction

    % Outlet
    u(end,2:end-1) = 2*u(end-1,2:end-1) - u(end-2,2:end-1);
    v(end,2:end-1) = 2*v(end-1,2:end-1) - v(end-2,2:end-1);
    P(end,2:end-1) = 2*P(end-1,2:end-1) - P(end-2,2:end-1);
    T(end,2:end-1) = 2*T(end-1,2:end-1) - T(end-2,2:end-1);

    % Wall BC will need serious consideration for pressure and temperature
    % variables, since they are extrapolated

    % Wall
    u(:,1) = 0; % Accounts for leading edge as well
    v(2:end,1) = 0;
    T(2:end,1) = T_inf;

    for i = 2:size(P,1)
        P(i,1) = 2*P(i,2) - P(i,3); % Wrong for now
    end
    
    % Apply changes to the primitive variables to the conservatives
    rho = P./(R*T);
    U = prim2cons(rho,u,v,T,cv);
end