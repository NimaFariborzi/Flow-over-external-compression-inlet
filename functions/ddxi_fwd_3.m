function dfdxi = ddxi_fwd_3(f,dxi,varargin)
%DDX_FWD_3 Computes first derivative of f in xi using first-order forward
%differencing for state-space arrays of 3 dimensions.
%   DDX_FWD_3(f,dxi) computes the partial derivative of 3D array f in xi with 
%   constant grid spacing dxi between points.

    % Initialize dfdxi array to same size as f
    dfdxi = zeros(size(f));
    % Compute dfdxi with first-order forward difference at all points up to the
    % right edge of the domain
    dfdxi(:,1:end-1,:) = (f(:,2:end,:) - f(:,1:end-1,:))/dxi;
    % Use backward difference at the right edge
    dfdxi(:,end,:) = (f(:,end,:) - f(:,end-1,:))/dxi;

    if nargin < 4
        return
    end
    cowl_rows = varargin{1};
    cowl_cols = varargin{2};
    % Use backward differencing on cowl left edge, if location provided
    i = cowl_rows(1);
    for j = cowl_cols
        dfdxi(:,i,j) = (f(:,i,j) - f(:,i-1,j))/dxi;
    end

end