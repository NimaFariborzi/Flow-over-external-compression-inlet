function dfdxi = ddxi_bwd(f,dxi,varargin)
%DDX_BWD Computes first derivative of f in xi using first-order backward
%differencing.
%   DDX_BWD(f,dxi) computes the partial derivative of 2D array f in xi with 
%   constant grid spacing dxi between points.

    % Initialize dfdxi array to same size as f
    dfdxi = zeros(size(f));
    % Compute dfdxi with first-order backward difference at all points to the
    % right of the left edge of the domain
    dfdxi(2:end,:) = (f(2:end,:) - f(1:end-1,:))/dxi;
    % Use forward difference at the left edge
    dfdxi(1,:) = (f(2,:) - f(1,:))/dxi;

    % Actually no difference due to cowl for ddxi_bwd since already using
    % backward differencing
end