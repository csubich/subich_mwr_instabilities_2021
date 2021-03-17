
function [depart_u, depart_phi] = semilag_traj(u,u_phi,x_phi,x_u,dt,semilag_fn,...
                                               in_depart_u,in_depart_phi)
% Calculate semilag departure points    
    if (nargin >= 7)
        depart_u = in_depart_u;
        depart_phi = in_depart_phi;
    else
        depart_u = x_u - dt*u; 
        depart_phi = x_phi - dt*u_phi;
    end
    resid_traj = 1.0;
    iter = 0;
    while resid_traj > 1e-8
        % Iterate via Trapezoidal rule until converged
        % Save old values
        depart_u_old = depart_u;
        depart_phi_old = depart_phi;

        % Perform trapezoidal iteration
        depart_u = x_u - 0.5*dt*(u + semilag_fn(depart_u,x_u)*u);
        depart_phi = x_phi - 0.5*dt*(u_phi + semilag_fn(depart_phi,x_u)*u);
        
        resid_traj = max(max(abs(depart_u-depart_u_old)),...
                         max(abs(depart_phi-depart_phi_old)));
        iter = iter + 1;
        if (iter > 10); break; end
    end
    return
end % End nested function