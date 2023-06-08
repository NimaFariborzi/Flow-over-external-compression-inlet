function d2fdet2 = d2det2(f,det)
%D2DY2 Computes second derivative of f in eta using second-order central
%differencing.
%   D2DY2(f,det) computes the second partial derivative of 2D array f in eta with 
%   constant grid spacing det between points.

    % Initialize d2fdet2 array to same size as f
    d2fdet2 = zeros(size(f));
    % Compute d2fdet2 with second-order central difference at all points in the
    % interior of the domain
    d2fdet2(:,2:end-1) = (f(:,3:end) - 2*f(:,2:end-1) + f(:,1:end-2))/(det^2);
    % Use second-order forward difference at the bottom edge
    d2fdet2(:,1) = (2*f(:,1) - 5*f(:,2) + 4*f(:,3) - f(:,4))/(det^2);
    % Use second-order backward difference at the top edge
    d2fdet2(:,end) = (2*f(:,end) - 5*f(:,end-1) + 4*f(:,end-2) - f(:,end-3))/(det^2);
end