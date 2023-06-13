function dfdet = ddet_bwd(f,det,varargin)
%DDY_BWD Computes first derivative of f in eta using first-order backward
%differencing.
%   DDY_BWD(f,det) computes the partial derivative of 2D array f in eta with 
%   constant grid spacing det between points.

    % Initialize dfdet array to same size as f
    dfdet = zeros(size(f));
    % Compute dfdet with first-order backward difference at all points above the
    % bottom edge of the domain
    dfdet(:,2:end) = (f(:,2:end) - f(:,1:end-1))/det;
    % Use forward difference at the bottom edge
    dfdet(:,1) = (f(:,2) - f(:,1))/det;

    if nargin < 4
        return
    end
    cowl_rows = varargin{1};
    cowl_cols = varargin{2};
    % Use one-sided differencing on cowl top, if location provided
    j2 = cowl_cols(2); % Top
    for i = cowl_rows
        dfdet(i,j2) = (f(i,j2+1) - f(i,j2))/det;
    end

end