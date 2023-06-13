function dfdxi = ddxi_central(f,dxi,varargin)
%DDX_CENTRAL Computes first derivative of f in xi using second-order central
%differencing.
%   DDX_CENTRAL(f,dxi) computes the partial derivative of 2D array f in xi with 
%   constant grid spacing dxi between points.

    % Initialize dfdxi array to same size as f
    dfdxi = zeros(size(f));
    % Compute dfdxi with second-order central difference at all points in the
    % interior of the domain
    dfdxi(2:end-1,:) = (f(3:end,:) - f(1:end-2,:))/(2*dxi);
    % Use second-order forward difference at the left edge
    dfdxi(1,:) = (-3*f(1,:) + 4*f(2,:) - f(3,:))/(2*dxi);
    % Use second-order backward difference at the right edge
    dfdxi(end,:) = (3*f(end,:) - 4*f(end-1,:) + f(end-2,:))/(2*dxi);

    if nargin < 4
        return
    end
    cowl_rows = varargin{1};
    cowl_cols = varargin{2};
    % Use backward differencing on cowl left edge, if location provided
    i = cowl_rows(1);
    for j = cowl_cols
        dfdxi(i,j) = (3*f(i,j) - 4*f(i-1,j) + f(i-2,j))/(2*dxi);
    end

end