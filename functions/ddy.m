function dfdy = ddy(f,xi_y,et_y,dxi,det,xi_dir,et_dir,varargin)
%DDY Computes first derivative of f in y on transformed grids.
%   DDY(f,xi_y,et_y,d_xi,d_et,xi_dir,et_dir) computes the partial derivative of
%   2D array f in y using grid transformations defined by xi_y and et_y with
%   grid spacing d_xi and d_et, taking derivatives in the direction specified by
%   the strings xi_dir and et_dir
%
%   xi_dir and et_dir must be one of the following: 'fwd', 'bwd', or 'cnt'

    % Select function for differencing in xi
    switch xi_dir
        case 'fwd'
            ddxi = @ddxi_fwd;
        case 'bwd'
            ddxi = @ddxi_bwd;
        case 'cnt'
            ddxi = @ddxi_central;
    end

    % Select function for differencing in eta
    switch et_dir
        case 'fwd'
            ddet = @ddet_fwd;
        case 'bwd'
            ddet = @ddet_bwd;
        case 'cnt'
            ddet = @ddet_central;
    end

    % Perform differencing
    dfdy = xi_y .* ddxi(f,dxi,varargin{:}) + et_y .* ddet(f,det,varargin{:});

end