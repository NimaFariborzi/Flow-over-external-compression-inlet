% Function to apply BC's to primitive variables
function [u, v, P, T, U] = apply_BCs_inlet(u, v, P, T, R, cv, u_inf, P_inf, T_inf, X, Y, cowl_rows, cowl_cols)

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
    
    % What's going on here is a bit convoluted but will be described in the
    % report
    for i = 2:size(P,1)
        % Get location of point ON wall
        x1 = X(i,1);
        y1 = Y(i,1);
        % Get location of point two ABOVE wall
        x4 = X(i,3);
        y4 = Y(i,3);
        % Get location of point two ABOVE wall and one LEFT of wall
        x3 = X(i-1,3);
        y3 = Y(i-1,3);
        % Choose appropriate y scale as distance to two points ABOVE wall
        ys = y4 - y1;
        % Get location on line normal from wall
        ang = atan2d(y1 - Y(i-1,1),x1 - X(i-1,1));
        x2 = x1 - ys*tand(ang);
        y2 = y1 + ys;
        % Find intersection point 2
        [ix2, iy2] = line_seg_intersect(x1,y1,x2,y2,x3,y3,x4,y4);
        % Interpolate pressure at intersection point 2
        d34 = norm([x4 - x3, y4 - y3]);
        d3i = norm([ix2 - x3, iy2 - y3]);
        P2 = P(i-1,3) + d3i/d34*(P(i,3) - P(i-1,3));
        % Find intersection point 1
        x4 = X(i,2);
        y4 = Y(i,2);
        x3 = X(i-1,2);
        y3 = Y(i-1,2);
        [ix1, iy1] = line_seg_intersect(x1,y1,x2,y2,x3,y3,x4,y4);
        % Interpolate pressure at intersection point 1
        d34 = norm([x4 - x3, y4 - y3]);
        d3i = norm([ix1 - x3, iy1 - y3]);
        P1 = P(i-1,2) + d3i/d34*(P(i,2) - P(i-1,2));
        % Get distance from wall to intersection point 1 and distance between
        % intersection points
        d1 = norm([ix1 - x1, iy1 - y1]);
        d2 = norm([ix1 - ix2, iy1 - iy2]);
        % Extrapolate pressure to wall
        P(i,1) = P1 - d1/d2*(P2 - P1);
    end

    % Apply BC's on cowl
    u(cowl_rows,cowl_cols) = 0;
    v(cowl_rows,cowl_cols) = 0;
    T(cowl_rows,cowl_cols) = 500; % Why?? Because it crashes if you try to
    % extrapolate or set adiabatic, and this is roughly the temp at which it is
    % adiabatic anyway...

    % Pressure on underside of cowl gets extrapolated below
    j1 = cowl_cols(1);
    for i = cowl_rows(2:end-1)
        d1 = norm([X(i,j1)-X(i,j1-1), Y(i,j1)-Y(i,j1-1)]);
        d2 = norm([X(i,j1-1)-X(i,j1-2), Y(i,j1-1)-Y(i,j1-2)]);
        P(i,j1) = P(i,j1-1) - d2/d1*(P(i,j1-2) - P(i,j1-1));
    end
    % Pressure on topside of cowl gets extrapolated above
    j2 = cowl_cols(2);
    for i = cowl_rows(2:end-1)
        d1 = norm([X(i,j2)-X(i,j2+1), Y(i,j2)-Y(i,j2+1)]);
        d2 = norm([X(i,j2+1)-X(i,j2+2), Y(i,j2+1)-Y(i,j2+2)]);
        P(i,j2) = P(i,j2+1) - d2/d1*(P(i,j2+2) - P(i,j2+1));
    end
    
    % Apply changes to the primitive variables to the conservatives
    rho = P./(R*T);
    U = prim2cons(rho,u,v,T,cv);
end