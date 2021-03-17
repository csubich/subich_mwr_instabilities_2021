function [L, lhs_phi, rhs_phi, lhs_u, rhs_u] = get_linearized_op(phi,u,hill_phi,x_phi,x_u,dt,...
                            imat_phi, imat_dphi, imat_u, imat_du, imat_u_phi,...
                            Dx_phi, Dx_u)
    % Get L, the linearized time-stepping operator corresponding to an implicit
    % trapezoidal discretization of the shallow water equations, linearized
    % about the presumed steady-state (phi,u)

    % Equation reminders:

    % Dphi/Dt + phi u_x = 0, or
    % phia + 0.5*dt*(phia ua_x) = ADV_phi(phib - 0.5*dt*(phib ub_x))
    
    % Du/Dt + g(H+phi)_x = 0, or
    % ua + 0.5*g*dt*(H+phia)_x = ADV_u(ub - 0.5*g*dt*(H+phib)_x)
    
    % Grid parameter:
    Nx = length(phi);
    global g; % Gravity
    
    % Build block matrices necessary to assemble the two-variable system
    % into a single matrix of size 2N x 2N
    I = eye(Nx);
    Z = zeros(Nx);
    Iphi = [I Z]; Iu = [Z I]; % Iphi * [phi; u] = [phi;0]
    
    % We need the derivative of the advective operator with respect to
    % perturbation u (u').  This depends on its value both before and after
    % the timestep, so it splits into before and after components
    
    % In general, dr/du | u=ubar is 
    % dr/du = dt/(2 + dt*ubar_x(x-rbar)) * (u'a + u'b(x - rbar)), so
    
    dubar_dx_u = imat_du*u;
    dubar_dx_phi = imat_dphi*(imat_u_phi*u);
    
    dra_u = (dt./(2 + dt*dubar_dx_u)).*Iu;
    drb_u = (dt./(2 + dt*dubar_dx_u)).*[Z imat_u];
    
    dra_phi = (dt./(2 + dt*dubar_dx_phi)).*[Z imat_u_phi];
    drb_phi = (dt./(2 + dt*dubar_dx_phi)).*[Z imat_phi*imat_u_phi];
    
    % The marginal impact of u' on ADV(f) is then
    % (-imat_d*f)*(dra + drb), which can be moved to the appropriate sides
    % of the equation
    
    dubardx = Dx_u*u;  % dubar/dx, via finite difference
    adv_phi_mean = phi.*(1 - dt/2*dubardx); % Mean-flow term being advected
    
    % Terms on the LHS of the phi' equation (involving phia and ua)
    lhs_phi = (1 + dt/2*dubardx).*Iphi + (0.5*dt*phi).*(Dx_u*Iu) - ...
            (-imat_dphi*adv_phi_mean).*dra_phi;
	% Terms on the RHS (involving phib and ub)
    rhs_phi = (-imat_dphi*adv_phi_mean).*drb_phi + ...
        imat_phi*((1 - dt/2*dubardx).*Iphi - (0.5*dt*phi).*(Dx_u*Iu));
    
    % Terns on the LHS of u'
    adv_u_mean = u - g*dt/2*Dx_phi*(phi + hill_phi);
    lhs_u = Iu + g*dt/2*Dx_phi*Iphi - (-imat_du*adv_u_mean).*dra_u;
    rhs_u = imat_u*(Iu - 0.5*g*dt*Dx_phi*Iphi) + (-imat_du*adv_u_mean).*drb_u;
    
    L = [lhs_phi;lhs_u] \ [rhs_phi;rhs_u];