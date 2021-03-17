% Worksheet for timestepping

% Perturbation based on eigenspectra (pick unstable mode)

if (~exist('pert_eps','var'))
    pert_eps = 1e-2; % Amplitude of perturbation relative to amp_hill
end


% Get scaling factors for the perturbation; normalize by the difference
% between the steady-state IC and uniform profiles (i.e. relative to the
% impact of the hill)
% norm0_phi = norm(phio - phi0)/sqrt(Nx);
% norm0_u = norm(uo-u0)/sqrt(Nx);

if (~exist('random_seed','var'))
    random_seed = 0;
end

if (~exist('INIT','var') || ~strcmp(INIT,'mode'))
    % Initialize with a random perturbation
    rng(20201218 + random_seed); % Apply a fixed seed
    pert_phi = randn(Nx,1);
    pert_phi = pert_eps*amp_hill*pert_phi;
    %pert_u = pert_eps*amp_hill/phi0*c0*randn(Nx,1);
    pert_u = 0*x_u;
else
    % Initialize with the precomputed most unstable mode
    pert_phi = real(cn_eigvec(1:Nx,maxidx));
    pert_u = real(cn_eigvec(Nx+(1:Nx),maxidx));
    rms_phi = norm(pert_phi)/Nx;
    % Normalize so phi perturbation has maximum amplitude of 1
    
    scale_factor = max(pert_phi)/amp_hill;
    
    pert_phi = pert_phi*pert_eps/scale_factor;
    pert_u = pert_u*pert_eps/scale_factor;
end
if (~exist('ploteach','var'))
    ploteach = 1;
end

if (~exist('FINTIME','var'))
    FINTIME = 4; % hours
end
% Solver type:
% 'newton' -- Netwon method, linearization at each iteration
% 'helmholtz' -- Approximate linearization based on constant phi at rest
if (~exist('SOLVER','var'))
    SOLVER = 'newton'; 
end
if (~strcmp(SOLVER,'newton') & ~exist('LINphi','var'))
    LINphi = 1.1*phi0; % Linearization level for helmholtz solver
end

if (~exist('iter_nonlin_limit','var'))
    if (strcmp(SOLVER,'newton'))
        iter_nonlin_limit = 10;
    else
        iter_nonlin_limit = 3;
    end
end

if (~exist('print_stats','var'))
    print_stats = 1;
end

% myphi = phio + pert_eps * pert_phi;
% myu = uo + pert_eps * pert_u;

% Simple initialization
% myphi = phi0 - hill_phi + pert_phi;
% myu = u0 + pert_u; %u0.*phi0./(phi0-hill_phi) + pert_u;

% Initialization about computed background state
myphi = phio + pert_phi;
myu = uo + pert_u;




semilag_fn = @interp_fourier;
% semilag_fn = @interp_linear;

% Get differential operators
[Dx_phi, Dx_u] = fd_ops(Nx,Lx);
imat_u_phi = semilag_fn(x_phi,x_u); % Average from u to phi points for semilag

% Build block matrices necessary to assemble the two-variable system
% into a single matrix of size 2N x 2N
I = eye(Nx);
Z = zeros(Nx);
Iphi = [I Z]; Iu = [Z I]; % Iphi * [phi; u] = [phi;0]


% Starting conditions for the timestep ("after" values from the "last"
% timestep)
ua = myu;
phia = myphi;

% Initial O(1) trajectory guess
depart_u = x_u - dt*ua;
depart_phi = x_phi - dt*imat_u_phi*ua;

amps = []; % Norm of disturbance, normalized by steady-state solution less phi0
residlog = [];

% dtlist = [0.6 1.4]*dt;
if (~exist('dtlist','var')) 
    dtlist = [dt];
end
nowtime = 0;

times = [];
phi_amp = [];
phi_max = [];
residlog = [];

step = 0;
%for step = 1:ceil(FINTIME*3600/dt)
while nowtime < FINTIME*3600
    step = step + 1;
    ub = ua;
    phib = phia;
    dtnow = dtlist(1 + mod(step,length(dtlist)));
%     dtnow = (0.8+0.4*rand(1))*dt;
%     dtnow = dt;

% Perform a single timestep

% Compute RHS on-grid
rhs_phi_grid = phib - dtnow/2*phib.*(Dx_u*ub);
rhs_u_grid = ub - g*dtnow/2*(Dx_phi*(phib + hill_phi));

resid_nonlin = 1;
iter_nonlin = 0;


while resid_nonlin > 1e-7 && iter_nonlin < iter_nonlin_limit
    
    
    
    % Trajectory iteration -- begin with the last best guess
    resid_traj = dx;
    iter_traj = 0;
    while resid_traj > 1e-8 && iter_traj < 10
        depart_u_old = depart_u;
        depart_phi_old = depart_phi;
        
        % Trapezoidal rule calculation
        depart_u = x_u - 0.5*dtnow*(ua + semilag_fn(depart_u_old,x_u)*ub);
        depart_phi = x_phi - 0.5*dtnow*(imat_u_phi*ua + semilag_fn(depart_phi_old,x_u)*ub);
        
        resid_traj = max(max(abs(depart_u - depart_u_old)), ...
                         max(abs(depart_phi - depart_phi_old)));
        iter_traj = iter_traj + 1;
    end
    
    % Advect RHS (old time level) to new grid position
    rhs_phi = semilag_fn(depart_phi,x_phi)*rhs_phi_grid;
    rhs_u = semilag_fn(depart_u,x_u)*rhs_u_grid;
    
    
    % LHS phi: phia + dt/2(phia * uax)
    % LHS u  : ua + gdt/2(phia + H)x
    
    lhs_phi = phia + dtnow/2*(phia.*(Dx_u*ua));
    lhs_u = ua + g*dtnow/2*(Dx_phi*(phia + hill_phi));
    
    resid_phi = rhs_phi - lhs_phi;
    resid_u = rhs_u - lhs_u;
    
    % Update based on Jacobian:
    if (strcmp(SOLVER,'newton'))
        Jacphi = Iphi + dtnow/2*((Dx_u*ua).*Iphi + phia.*(Dx_u*Iu));
        Jacu   = Iu + g*dtnow/2*(Dx_phi*Iphi);
        prime = [Jacphi; Jacu] \ [resid_phi; resid_u];
        phiprime = prime(1:Nx);
        uprime = prime(Nx+(1:Nx));
    else % Helmholtz solver
        % linearized: phi' + dt/2(linphi*u'_x) = resid_phi
        %             u'  + gdt/2(phi'_x)    = resid_u
        %          Du u' + gdt/2 Du Dphi phi'= Du resid_u
        %          phi' - g*linphi dt^2 / 4 Du Dphi phi' = resid_phi - linphi dt/2*
        %                                                        Du resid_u
        Helm = I - g*LINphi*dtnow^2/4 * (Dx_u * Dx_phi);
        rhshelm = resid_phi - LINphi*dtnow/2 * (Dx_u*resid_u);
        phiprime = Helm \ rhshelm;
        uprime = resid_u - g*dtnow/2*(Dx_phi*phiprime);        
    end
    
    
    %resid_nonlin = max(abs(prime));
    resid_nonlin = Nx^-0.5*(norm(phiprime)/phi0 + norm(uprime)/u0);
    if (print_stats)
        fprintf('%d: %g\n',iter_nonlin+1, resid_nonlin);
    end
    
    phia = phia + phiprime;
    ua = ua + uprime;
    
    % Restore mass
%     phia = phia + mean(myphi) - mean(phia);
    
    iter_nonlin = iter_nonlin+1;
  
end


nowtime = nowtime + dtnow;
times(step) = nowtime;

%amps(step) = Nx^-0.5*(norm(phia-phio)/phi0 + norm(ua-uo)/u0);
phi_amp(step) = sqrt(mean((phia-phio).^2));
phi_max(step) = max(abs(phia-phio));
residlog(step) = resid_nonlin;

if (mod(step,ploteach) == 0)
    %plot(x_phi,phia-phi0,x_phi,phio-phi0);
    %plot(x_phi,phia-phio,x_phi,pert_phi);
    %plot(x_phi,phia+hill_phi,x_phi,phio+hill_phi);
    plot(x_phi,phia+hill_phi,x_phi,myphi+hill_phi);
    
    title(nowtime)
    drawnow; pause(0.01);
end


end