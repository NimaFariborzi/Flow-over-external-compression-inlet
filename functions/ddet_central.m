function dfdet = ddet_central(f,det)
%DDET_CENTRAL Computes first derivative of f in eta using second-order central
%differencing.
%   DDET_CENTRAL(f,det) computes the partial derivative of 2D array f in eta with 
%   constant grid spacing det between points.

    % Initialize dfdet array to same size as f
    dfdet = zeros(size(f));
    % Compute dfdet with second-order central difference at all points in the
    % interior of the domain
    dfdet(:,2:end-1) = (f(:,3:end) - f(:,1:end-2))/(2*det);
    % Use second-order forward difference at the bottom edge
    dfdet(:,1) = (-3*f(:,1) + 4*f(:,2) - f(:,3))/(2*det);
    % Use second-order backward difference at the top edge
    dfdet(:,end) = (3*f(:,end) - 4*f(:,end-1) + f(:,end-2))/(2*det);
end