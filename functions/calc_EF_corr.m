% Function to calculate the values of E and F to be used during corrector step
% calculations
function [E, F] = calc_EF_corr(U,u,v,P,T,cp,Pr,dx,dy)
    
    % NEEDS MODIFICATION

    % Allocate E and F
    E = zeros(size(U));
    F = zeros(size(U));
    % Calculate viscosity and thermal conductivity
    mu = sutherland(T);
    k = mu*cp/Pr;
    % Calculate normal stresses
    tau_xx = mu.*(4/3*ddxi_fwd(u,dx) - 2/3*ddet_central(v,dy));
    tau_yy = mu.*(4/3*ddet_fwd(v,dy) - 2/3*ddxi_central(u,dx));
    % Calculate shear stresses (different for E and F due to finite differences)
    tau_xy_E = mu.*(ddet_central(u,dy) + ddxi_fwd(v,dx));
    tau_xy_F = mu.*(ddxi_central(v,dx) + ddet_fwd(u,dy));
    % Calculate heat fluxes
    q_x = -k.*ddxi_fwd(T,dx);
    q_y = -k.*ddet_fwd(T,dy);
    % First slice
    E(1,:,:) = U(2,:,:);
    F(1,:,:) = U(3,:,:);
    % Second slice
    E(2,:,:) = squeeze(U(2,:,:)).*u + P - tau_xx;
    F(2,:,:) = squeeze(U(2,:,:)).*v     - tau_xy_F;
    % Third slice
    E(3,:,:) = squeeze(U(3,:,:)).*u     - tau_xy_E;
    F(3,:,:) = squeeze(U(3,:,:)).*v + P - tau_yy;
    % Fourth slice
    E(4,:,:) = (squeeze(U(4,:,:)) + P).*u - u.*tau_xx - v.*tau_xy_E + q_x;
    F(4,:,:) = (squeeze(U(4,:,:)) + P).*v - v.*tau_yy - u.*tau_xy_F + q_y;
end