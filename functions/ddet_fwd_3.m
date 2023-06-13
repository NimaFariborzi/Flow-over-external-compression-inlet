function dfdet = ddet_fwd_3(f,det,varargin)
%DDET_FWD_3 Computes first derivative of f in eta using first-order forward
%differencing for state-space arrays of 3 dimensions.
%   DDET_FWD_3(f,det) computes the partial derivative of 3D array f in eta with 
%   constant grid spacing det between points.

    % Initialize dfdet array to same size as f
    dfdet = zeros(size(f));
    % Compute dfdet with first-order forward difference at all points up to the
    % top edge of the domain
    dfdet(:,:,1:end-1) = (f(:,:,2:end) - f(:,:,1:end-1))/det;
    % Use backward difference at the top edge
    dfdet(:,:,end) = (f(:,:,end) - f(:,:,end-1))/det;

    if nargin < 4
        return
    end
    cowl_rows = varargin{1};
    cowl_cols = varargin{2};
    % Use one-sided differencing on cowl bottom, if location provided
    j1 = cowl_cols(1); % Bottom
    for i = cowl_rows
        dfdet(:,i,j1) = (f(:,i,j1) - f(:,i,j1-1))/det;
    end

end