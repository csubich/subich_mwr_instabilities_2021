function [imat, imat_d] = interp_fourier(xi,x_grid)
    % Build interpolating matrix to points xi, from a given grid
    % Ref: Boyd, appendix F2
    Nx = length(x_grid);
    Lx = (max(x_grid) - min(x_grid))*(Nx)/(Nx-1);
    xf = (xi - x_grid(1))/Lx*2*pi;
    xf = xf(:);
    if (mod(Nx,2) == 0) % Even case
        Nhalf = Nx/2;
        xc = 2*pi*(0:Nx-1)/Nx;
        imat = 1/(2*Nhalf).*sin(Nhalf*(xf-xc)).*cot(0.5*(xf-xc));
        if (nargout > 1)
            imat_d = 1/(2*Nhalf).*(Nhalf.*cos(Nhalf*(xf-xc)).*cot(0.5*(xf-xc)) + ...
                                   -0.5*sin(Nhalf*(xf-xc))./sin(0.5*(xf-xc)).^2);
        end
    else % Odd case -- work on a doubled grid
        Nhalf = Nx;
        xc = 2*pi*(0:(2*Nx-1))/(2*Nx);
        % Cardinal functions are the sum of the left and right halves
        imat = 1/(2*Nhalf).*sin(Nhalf*(xf/2-xc(1:Nx))).*cot(0.5*(xf/2-xc(1:Nx))) + ...
           1/(2*Nhalf).*sin(Nhalf*(xf/2-xc(Nx+(1:Nx)))).*cot(0.5*(xf/2-xc(Nx+(1:Nx))));

        if (nargout > 1)
            imat_d = 1/(2*Nhalf).*(Nhalf/2*cos(Nhalf*(xf/2-xc(1:Nx))).*cot(0.5*(xf/2-xc(1:Nx))) + ...
                                   -0.25*sin(Nhalf*(xf/2-xc(1:Nx)))./sin(0.5*(xf/2-xc(1:Nx))).^2) + ...
                     1/(2*Nhalf).*(Nhalf/2*cos(Nhalf*(xf/2-xc(Nx+(1:Nx)))).*cot(0.5*(xf/2-xc(Nx+(1:Nx)))) + ...
                                   -0.25*sin(Nhalf*(xf/2-xc(Nx+(1:Nx))))./sin(0.5*(xf/2-xc(Nx+(1:Nx)))).^2);
        end
    end
    
    % Replace any nans with 1, as the limit of sin*cot
    imat(isnan(imat)) = 1;
    
    % Set row sum of interpolating matrix to 1
    imat = imat - diag(sum(imat,2)-1);
    if (nargout > 1)
        imat_d(isnan(imat_d)) = 0;
        imat_d = imat_d*2*pi/Lx; % Scale by inverse grid length
        % Set row sum of differentation matrix to 0
        imat_d = imat_d - diag(sum(imat_d,2));
    end
    
end
