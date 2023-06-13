function dfdx = ddx(f,xi_x,et_x,dxi,det,xi_dir,et_dir,varargin)
%DDX Computes first derivative of f in x on transformed grids.
%   DDX(f,xi_x,et_x,d_xi,d_et,xi_dir,et_dir) computes the partial derivative of
%   2D array f in x using grid transformations defined by xi_x and et_x with
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
    dfdx = xi_x .* ddxi(f,dxi,varargin{:}) + et_x .* ddet(f,det,varargin{:});

end