% Test script for initial work re: off-centering

% This script:
% * builds a linearized initial condition for a single physical
% configuration
% * builds the trapezoidal operators corresponding to this initial
% condition
% * plots the resulting spectrum

%% physical paramters
%phi0 = 10; g = 10; u0 = 1;
if (~exist('phi0','var'))
    phi0 = 6000; u0 = 20;
    global g;
    g = 9.81;
    c = sqrt(phi0*g); % long wave phase speed



    %% numerical parameters
    %dx = 1; Nx = 129; Lx = dx*Nx;
    dx = 1000; Nx = 257; Lx = dx*Nx;
    [Dx_phi, Dx_u] = fd_ops(Nx,Lx); % Get finite-difference operators

    % Define the grid:
    x_phi = Lx*(-0.5 + (-0.5 + (1:Nx)')/Nx);
    x_u = x_phi + dx/2;

    % Temporal parameters
    %cfl = 1.2; 
    %cfl = 1.18;
    cfl = 1.683;
    dt = cfl*dx/u0;
    
    %% hill parameters (defined in terms of dx units)
    %hill_width = 5*dx; % peak-to-peak wavelength of topography
    %hill_width = 257e3/39;
    hill_width = 257e3/50;
    % Limit k_hill so that a fixed number of periods fit in the domain
    num_hill = floor(Lx/hill_width);
    k_hill = 2*pi*num_hill/Lx;
end
%% Off-centering parameter
    if (~exist('alpha','var'))
        alpha = 0.50; % Off-centering parameter
end

%amp_hill = 1e-2*phi0;
if (~exist('amp_hill','var'))
    amp_hill = 1.667e-2*phi0;
end

hill_phi = amp_hill*sin(k_hill*x_phi);

%% Get leading-order steady state

% The hill itself is defined as amp_hill*sin(k_hill*x), and we need to take
% it one Fourier mode at a time.

% Positive frequency:
slhat = exp(-1i*k_hill*u0*dt);
dxhat = 2i*dt/dx*sin(k_hill*dx/2);
hhat = amp_hill/2i; % (sin k = (exp(ik) - exp(-ik)) / 2i

phiuhat_pos = [1 - slhat, phi0*dxhat*(alpha + (1-alpha)*slhat) ;
           g*dxhat*(alpha + (1-alpha)*slhat), 1 - slhat] \ ...
           [0; -dxhat*g*(alpha + (1-alpha)*slhat)*hhat];

% Negative frequency

slhat = exp(1i*k_hill*u0*dt);
dxhat = -2i*dt/dx*sin(k_hill*dx/2);
hhat = -amp_hill/2i; % (sin k = (exp(ik) - exp(-ik)) / 2i

phiuhat_neg = [1 - slhat, phi0*dxhat*(alpha + (1-alpha)*slhat) ;
           g*dxhat*(alpha + (1-alpha)*slhat), 1 - slhat] \ ...
           [0; -dxhat*g*(alpha + (1-alpha)*slhat)*hhat];

% Add the basic, positive, and negative components to get the background
% state
phibar = phi0 + phiuhat_pos(1)*exp(1i*k_hill*x_phi) + phiuhat_neg(1)*exp(-1i*k_hill*x_phi);
ubar = u0 + phiuhat_pos(2)*exp(1i*k_hill*x_u) + phiuhat_neg(2)*exp(-1i*k_hill*x_u);

%% With phibar/ubar, define the semi-Lagrangian operator.

% Begin with displacements.  Since we have an analytic expression for ubar
% above, we can use that for delta-r displacements as well.  This
% implicitly assumes the use of Fourier interpolation in the
% semi-Lagrangian operator, where the +-k_hill frequency waves are
% interpolated exactly without amplitude changes.

delr0 = u0*dt; % Displacement from the flat-topography state (u0)
% Additional displacements due to ubar (exclusive of u0)
% u-points
delrbar_u = dt*(alpha*(ubar-u0) + ... % Arrival point
           (1-alpha)*(phiuhat_pos(2)*exp(1i*k_hill*(x_u - delr0)) + ...
                      phiuhat_neg(2)*exp(-1i*k_hill*(x_u-delr0)))); % Departure point
% phi-points
delrbar_phi = dt*(alpha*(phiuhat_pos(2)*exp(1i*k_hill*x_phi) + ...
                         phiuhat_neg(2)*exp(-1i*k_hill*x_phi)) + ...
                  (1-alpha)*(phiuhat_pos(2)*exp(1i*k_hill*(x_phi-delr0)) + ...
                         phiuhat_neg(2)*exp(-1i*k_hill*(x_phi-delr0))));
                       
% Semi-Lagrangian interpolation operator and the implied derivative, which
% is based on delr0
% [sl0_u, dsl0_u]     = interp_fourier(x_u-delr0, x_u);
% [sl0_phi, dsl0_phi] = interp_fourier(x_phi-delr0, x_phi);
% Interpolation from u to phi grid (for ua' terms in phi equations)
[imat_u_phi] = interp_fourier(x_phi,x_u);

%% Build the linearized operator, largely following get_linearized_op 

% Building-block matrices
I = eye(Nx); Z = zeros(Nx);
Iphi = [I Z]; Iu = [Z I];

% Incremental semi-lagrangian operators

% semilag-bar, from the mean flow
% slbar_u = -(delrbar_u.*dsl0_u); % (delrbar * f_x(x-delr0)
% slbar_phi = -(delrbar_phi.*dsl0_phi);

[sl0_u, dsl0_u] = interp_fourier(x_u - (delr0 + delrbar_u), x_u);
[sl0_phi, dsl0_phi] = interp_fourier(x_phi - (delr0 + delrbar_phi), x_phi);
slbar_u = 0; slbar_phi = 0;

% In a compact matrix contation, semilag-prime needs to split: 
% slp(f) = -fx(x-delr0)*u' = (-dsl0*f).*{alpha, 1-alpha}Iu
phi_lhs = Iphi + alpha*dt*(phibar.*(Dx_u*Iu) + ...
                          (Dx_u*ubar).*Iphi) + ...
          (dsl0_phi*(phibar - (1 - alpha)*dt*phi0*Dx_u*ubar)).*(alpha*dt*(imat_u_phi*Iu));
% phi_rhs = (sl0_phi + slbar_phi)*(Iphi - (1-alpha)*dt*phi0*Dx_u*Iu) + ...
%            sl0_phi*(-(1-alpha)*dt)*((phibar-phi0).*(Dx_u*Iu)) + ...
%           (dsl0_phi*(phibar - (1 - alpha)*phi0*Dx_u*ubar)).*((1-alpha)*(sl0_u*(imat_u_phi*Iu)));
phi_rhs = (sl0_phi + slbar_phi)*(Iphi - (1-alpha)*dt*phibar.*Dx_u*Iu) - ...
          (dsl0_phi*(phibar - (1 - alpha)*dt*phi0*Dx_u*ubar)).*((1-alpha)*dt*(sl0_phi*(imat_u_phi*Iu)));
       
u_lhs = Iu + alpha*g*dt*Dx_phi*Iphi + ...
        (dsl0_u*(ubar - (1-alpha)*g*dt*Dx_phi*(phibar + hill_phi))).*(alpha*dt*Iu);
u_rhs = (sl0_u + slbar_u)*(Iu - (1-alpha)*g*dt*Dx_phi*Iphi) - ...
        (dsl0_u*(ubar - (1-alpha)*g*dt*Dx_phi*(phibar + hill_phi))).*((1-alpha)*dt*sl0_u*Iu);
    
cn_me = [phi_lhs; u_lhs] \ [phi_rhs; u_rhs];
imat_u = sl0_u + slbar_u;

% [imat_u, imat_du] = interp_fourier(x_u - (delr0 + delrbar_u), x_u);
% [imat_phi, imat_dphi] = interp_fourier(x_phi - (delr0 + delrbar_phi), x_phi);
% imat_u = sl0_u + slbar_u;
% imat_du = dsl0_u;
% imat_phi = sl0_phi + slbar_phi;
% imat_dphi = dsl0_phi;
% [cn_lin, lhs_phi_lin, rhs_phi_lin, lhs_u_lin, rhs_u_lin] = ...
%     get_linearized_op(phibar,ubar,hill_phi,x_phi,x_u,dt,...
%                       imat_phi, imat_dphi, ...
%                       imat_u, imat_du, ...
%                       imat_u_phi, Dx_phi, Dx_u);