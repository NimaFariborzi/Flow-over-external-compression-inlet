% Function to calculate the values of E and F to be used during corrector step
% calculations
function [E, F] = calc_EF_corr(U,u,v,P,T,cp,Pr,d_xi,d_et,xi_x,et_x,et_y,xi_y,varargin)
    % Allocate E and F
    E = zeros(size(U));
    F = zeros(size(U));
    % Calculate viscosity and thermal conductivity
    mu = sutherland(T);
    k = mu*cp/Pr;
    % Calculate normal stresses
    tau_xx = mu.*(4/3*ddx(u,xi_x,et_x,d_xi,d_et,'fwd','cnt',varargin{:}) - 2/3*ddy(v,xi_y,et_y,d_xi,d_et,'fwd','cnt',varargin{:}));
    tau_yy = mu.*(4/3*ddy(v,xi_y,et_y,d_xi,d_et,'cnt','fwd',varargin{:}) - 2/3*ddx(u,xi_x,et_x,d_xi,d_et,'cnt','fwd',varargin{:}));
    % Calculate shear stresses (different for E and F due to finite differences)
    tau_xy_E = mu.*(ddy(u,xi_y,et_y,d_xi,d_et,'fwd','cnt',varargin{:}) + ddx(v,xi_x,et_x,d_xi,d_et,'fwd','cnt',varargin{:}));
    tau_xy_F = mu.*(ddx(v,xi_x,et_x,d_xi,d_et,'cnt','fwd',varargin{:}) + ddy(u,xi_y,et_y,d_xi,d_et,'cnt','fwd',varargin{:}));
    % Calculate heat fluxes
    q_x = -k.*ddx(T,xi_x,et_x,d_xi,d_et,'fwd','cnt',varargin{:});
    q_y = -k.*ddy(T,xi_y,et_y,d_xi,d_et,'cnt','fwd',varargin{:});
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