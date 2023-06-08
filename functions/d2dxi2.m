function d2fdxi2 = d2dxi2(f,dxi)
%D2DX2 Computes second derivative of f in xi using second-order central
%differencing.
%   D2DX2(f,dxi) computes the second partial derivative of 2D array f in xi with 
%   constant grid spacing dxi between points.

    % Initialize d2fdxi2 array to same size as f
    d2fdxi2 = zeros(size(f));
    % Compute d2fdxi2 with second-order central difference at all points in the
    % interior of the domain
    d2fdxi2(2:end-1,:) = (f(3:end,:) - 2*f(2:end-1,:) + f(1:end-2,:))/(dxi^2);
    % Use second-order forward difference at the left edge
    d2fdxi2(1,:) = (2*f(1,:) - 5*f(2,:) + 4*f(3,:) - f(4,:))/(dxi^2);
    % Use second-order backward difference at the right edge
    d2fdxi2(end,:) = (2*f(end,:) - 5*f(end-1,:) + 4*f(end-2,:) - f(end-3,:))/(dxi^2);
end