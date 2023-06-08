function dfdxi = ddxi_fwd(f,dxi)
%DDX_FWD Computes first derivative of f in xi using first-order forward
%differencing.
%   DDX_FWD(f,dxi) computes the partial derivative of 2D array f in xi with 
%   constant grid spacing dxi between points.

    % Initialize dfdxi array to same size as f
    dfdxi = zeros(size(f));
    % Compute dfdxi with first-order forward difference at all points up to the
    % right edge of the domain
    dfdxi(1:end-1,:) = (f(2:end,:) - f(1:end-1,:))/dxi;
    % Use backward difference at the right edge
    dfdxi(end,:) = (f(end,:) - f(end-1,:))/dxi;
end