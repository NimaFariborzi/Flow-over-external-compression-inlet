function dfdet = ddet_fwd(f,det)
%DDET_FWD Computes first derivative of f in eta using first-order forward
%differencing.
%   DDET_FWD(f,det) computes the partial derivative of 2D array f in eta with 
%   constant grid spacing det between points.

    % Initialize dfdet array to same size as f
    dfdet = zeros(size(f));
    % Compute dfdet with first-order forward difference at all points up to the
    % top edge of the domain
    dfdet(:,1:end-1) = (f(:,2:end) - f(:,1:end-1))/det;
    % Use backward difference at the top edge
    dfdet(:,end) = (f(:,end) - f(:,end-1))/det;
end