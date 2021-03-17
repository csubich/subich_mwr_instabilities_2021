function [imat, imat_d] = interp_cubic(xi, x_grid);
    % Build interpolating matrix to points xi, from a given uniform grid
    % x_grid
    
    Nx = length(x_grid);
    xmin = min(x_grid);
    Lx = (max(x_grid) - xmin)*Nx/(Nx-1);
    dx = Lx/Nx;
    
    Npt = length(xi);
    
    % Locate points xi in grid-index coordinates (base 0)
    grid_locations = mod((xi(:)-xmin)/dx,Nx);
    integer_locations = mod(floor(grid_locations),Nx);
    frac = grid_locations - integer_locations;
    
    imat = zeros(Npt,Nx);
    
    % Compute the 4-point Lagrangian interpolation coefficients
   
    % At integer_location - 1
    % (x-x1)(x-x2)(x-x3)/(x0-x1)(x0-x2)(x0-x3)
    pm1 = (frac-0).*(frac-1).*(frac-2)./((-1)*(-2)*(-3));
    
    % At integer_location
    % (x-x0)(x-x2)(x-x3)/(x1-x0)(x1-x2)(x1-x3)
    p0  = (frac+1).*(frac-1).*(frac-2)./(1*(-1)*(-2));
    
    % At integer_location+1
    % (x-x0)(x-x1)(x-x3)/(x2-x0)(x2-x1)(x2-x3)
    pp1 = (frac+1).*(frac-0).*(frac-2)./(2*1*(-1));
    
    % At integer_location+4
    % (x-x0)(x-x1)(x-x2)/(x3-x0)(x3-x1)(x3-x2)
    pp2 = (frac+1).*(frac-0).*(frac-1)./(3*2*1);
    
    imat(sub2ind([Npt,Nx],(1:Npt)',1+mod(-1+integer_locations,Nx))) = pm1;
    imat(sub2ind([Npt,Nx],(1:Npt)',1+mod(0 +integer_locations,Nx))) = p0;
    imat(sub2ind([Npt,Nx],(1:Npt)',1+mod( 1+integer_locations,Nx))) = pp1;
    imat(sub2ind([Npt,Nx],(1:Npt)',1+mod( 2+integer_locations,Nx))) = pp2;
    
    if (nargout > 1)
        imat_d = zeros(Npt,Nx);
        % Derivatives of the point formulas
        dm1 = ((frac-1).*(frac-2) + (frac-0).*(frac-2) + (frac-0).*(frac-1))/(-6*dx);
        d0  = ((frac-1).*(frac-2) + (frac+1).*(frac-2) + (frac+1).*(frac-1))/(2*dx);
        dp1 = ((frac-0).*(frac-2) + (frac+1).*(frac-2) + (frac+1).*(frac-0))/(-2*dx);
        dp2 = ((frac-0).*(frac-1) + (frac+1).*(frac-1) + (frac+1).*(frac-0))/(6*dx);
        
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+mod(-1+integer_locations,Nx))) = dm1;
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+mod(0 +integer_locations,Nx))) = d0;
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+mod( 1+integer_locations,Nx))) = dp1;
        imat_d(sub2ind([Npt,Nx],(1:Npt)',1+mod( 2+integer_locations,Nx))) = dp2;       
    end
        
    
end