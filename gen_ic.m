function [phi, u, resid, imat_u, imat_du, imat_phi, imat_dphi] = gen_ic(phib,ub,x_phi,x_u,...
                                   dt,semilag_fn,hill_phi,...
                                   Dx_phi,Dx_u,imat_u_phi,itercount,tol,regularize)
% GEN_IC -- generate steady-state initial conditions based on an idealized
% configuration
%  [phi, u, resid] = gen_ic(phib,ub,x_phi,x_u,dt,semilag_fn,hill_fn,Dx_phi,Dx_u,itercount)
%  iteratively optimizes the (phib, ub) fields such that the Crank-Nicholson
%  problem has a steady-state solution of [phi, u].  For full generality,
%  the semilag operator is passed as a function that is expected to return
%  [imat, imat_d], where (imat) is an interpolation matrix and (imat_d) is
%  the implied spatial derivative (i.e. imat*(d/dx f) = imat_d(f), noting
%  that the derivative is taken in the interpolation function space.)

% Grid size
Nx = length(x_u);
global g

% Preparatory work: find u at the phi points for trajectory calculations
u = ub;
phi = phib;
if (~exist('tol'))
    tol = 1e-4;
end
if (~exist('regularize'))
    regularize = 0;
end
% tol

% Matrix for interpolating from u-grid to phi-grid
%[imat_u_phi] = semilag_fn(x_phi, x_u);

% Necessary helper matrices for building the Jacobian
% Full vector will be stacked as [phi;u]
I = speye(Nx);
Z = sparse(Nx,Nx);
Iphi = [I Z]; Iu = [Z I];

% Mean values of u and phi
mean_u = mean(u);
mean_phi = mean(phi);

for iter=0:itercount
    % Get trajectories
    u_phi = imat_u_phi*u;
    if (iter == 0)
        [depart_u, depart_phi] = semilag_traj(u,u_phi,x_phi,x_u,dt,semilag_fn);
    else
        % Re-use departure points as initial guess
        [depart_u, depart_phi] = semilag_traj(u,u_phi,x_phi,x_u,dt,semilag_fn,...
                            depart_u,depart_phi);
    end
    
    % Build semilag operators
    [imat_u, imat_du] = semilag_fn(depart_u,x_u);
    [imat_phi, imat_dphi] = semilag_fn(depart_phi,x_phi);
    
    % Get residuals in phi and u equations
    [resid_phi, resid_u] = get_resid(phi,u,hill_phi,g,dt,Dx_phi,Dx_u,...
                                     imat_phi, imat_u);
                                 
% 	fprintf('%7.5g %7.5g \n',norm(resid_phi),norm(resid_u))
    if (norm(resid_phi) + norm(resid_u) < tol)
        break
    end
                                 
	if (iter ~= itercount)
        % Build Jacobian
        
        % Note that the semilag matrix has a dependency on the trajectory
        % length and thus u.  Trajectories are given by r = -(u +
        % imat*u)*0.5*dt, so dr = -(du + imat*du + dr*ux) (plus h.o.t.), or
        dr_u = -dt/2*diag(1.0./(1+dt/2*imat_du*u))*(I+imat_u);
        %u_phi = imat_u_phi*u;
        dr_phi = -dt/2*diag(1.0./(1+dt/2*imat_dphi*u_phi))*(imat_u_phi + imat_phi);
        
        % Thus, the derivative of imat*f is imat*df + (imat_d*f)*(dr*du)
        
        % Associated variables
        hphi = hill_phi + phi;
        dxu = Dx_u*u; % Lives at phi points
        dxhphi = Dx_phi*hphi; % Lives at u points
        
        % Phi residual:
        % resid_phi = (phi - imat_phi*phi)/dt + 0.5*(imat_phi*(phi.*dxu) + (phi.*dxu));
        
        DrDphi = (Iphi - (imat_phi*Iphi + (imat_dphi*phi).*(dr_phi*Iu)))/dt + ...
                 0.5*(phi.*(Dx_u*Iu) + dxu.*Iphi + ...
                      imat_phi*(phi.*(Dx_u*Iu) + dxu.*Iphi) + ...
                      (imat_dphi*(phi.*dxu)).*(dr_phi*Iu));
                  
        % U residual:
        % resid_u = (u - imat_u*u)/dt + 0.5*g*(dxhphi + imat_u*dxhphi);
        
        DrDu = (Iu - (imat_u*Iu + (imat_du*u).*(dr_u*Iu)))/dt + ...
                0.5*g*(Dx_phi*Iphi + imat_u*(Dx_phi*Iphi) + ...
                    (imat_du*dxhphi).*(dr_u*Iu));
        
        % The Jacobian as written has two zero eigenvalues by analogy to f_t =
        % f_x with periodic boundary conditions, so augment the system s.t. we
        % preserve the mean of phi and u
        
        Jac = [DrDphi ones(Nx,1) zeros(Nx,1);
               DrDu   zeros(Nx,1) ones(Nx,1);
               ones(1,Nx)/Nx zeros(1,Nx+2);
               zeros(1,Nx) ones(1,Nx)/Nx 0 0];
               %zeros(1,Nx) 1/Nx*1./u' 0 0];
           
%         fprintf(' %7.5g', cond(Jac))
        
        % Get correction
        %corr = Jac \ [-resid_phi; -resid_u; 0;mean(1./u)-mean(1./ub)];
        corr = Jac \ [-resid_phi; -resid_u; 0;0];
        if (regularize)
            %steep = -[resid_phi' resid_u' 0 mean(1./u)-mean(1./ub)]*Jac; % Get direction of steepest descent
            steep = -[resid_phi' resid_u' 0 0]*Jac; % Get direction of steepest descent
            % Limit the step by the cosine similarity between the two
            relax = steep*corr/norm(steep)/norm(corr);
        else
            relax = 1;
        end
        phi = phi + relax*corr(1:Nx);
        u = u + relax*corr(Nx+(1:Nx));
        
%         fprintf(' %.5e\n',relax);
        % Adjust phi, u
    else
%         fprintf('\n')
    end

end

resid = norm(resid_u) + norm(resid_phi);

end % End function

function [resid_phi, resid_u] = get_resid(phi, u, hill_phi, g, dt, Dx_phi, Dx_u, imat_phi, imat_u)
    % Shallow-water equations are:
    
    % dphi/dt + (phi u)_x = 0
    % du/dt + (u^2/2 + g(H+phi))_x = 0
    
    % In advective form, this is:
    % Dphi/Dt + phi u_x = 0
    % Du/Dt + g(H+phi)_x = 0
    
    % Get augmented H+phi variable
    hphi = hill_phi + phi;
    
    % Note that:
    % Df/Dt -> (f - f_d)/dt
    % (f)_x -> 0.5*(f_x + (f_x)_d)
    
    % That is, we evaluate linear terms before advecting. 
    
    dxu = Dx_u*u; % Lives at phi points
    dxhphi = Dx_phi*hphi; % Lives at u points
    
    % Evaluate residuals
    resid_phi = (phi - imat_phi*phi)/dt + 0.5*(imat_phi*(phi.*dxu) + (phi.*dxu));
    resid_u = (u - imat_u*u)/dt + 0.5*g*(dxhphi + imat_u*dxhphi);
    
end % End nested function
    
