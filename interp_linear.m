function [imat, imat_d] = interp_linear(xi, x_grid)
    % Build interpolating matrix to points xi, from a given grid
    Nx = length(x_grid);
    xmin = min(x_grid);
    Lx = (max(x_grid) - min(x_grid))*Nx/(Nx-1);
    dx = Lx/Nx;
    
    Npt = length(xi);
    
    % Compute locations xi in grid coordinates (base 0)
    grid_locations = mod((xi-xmin)/dx,Nx);
    integer_locations = mod(floor(grid_locations),Nx);
    fraction_locations = grid_locations - integer_locations;
    
    imat = zeros(Npt,Nx);
    
    % In parallel, compute the 2-point Lagrangian interpolation
    % coefficients
    lhs_coef = 1-fraction_locations;
    rhs_coef = fraction_locations;
    
    imat(sub2ind([Npt,Nx],(1:Npt)',1+integer_locations)) = lhs_coef;
    imat(sub2ind([Npt,Nx],(1:Npt)',1+mod(1+integer_locations,Nx))) = rhs_coef;
    
    if (nargout > 1)
        imat_d = zeros(Npt,Nx);
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+integer_locations)) = 1/dx;
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+mod(1+integer_locations,Nx))) = -1/dx;
    end
end