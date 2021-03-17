% Script to generate instability curves for the paper

dx = 1000;
Nx = 513;
% Nx = 257;
Lx = dx*Nx;
u0 = 20; % Mean velocity
phi0 = 6000; % Mean depth
global g
g = 9.81;
c = sqrt(g*phi0); % phase speed

% Hill amplitude
% 1.67e-2 -> 100m hills on 6km phi0
amps = 1.667e-2;%logspace(-4,-2,15);1e-3; 

% Hill parameters
periodlist = 5:floor(Nx/3); % Number of periods of topography in domain

% CFL parameters
Ncfl = 200;
% Ncfl = 32;
cflmin = 0.05; cflmax = 4;
cfl_list = cflmin + (cflmax-cflmin)*linspace(0,1,Ncfl);
Ncfl = length(cfl_list);


maxcn_save = zeros(length(cfl_list),length(periodlist),length(amps));
maxadv_save = zeros(length(cfl_list),length(periodlist),length(amps));
ic_resids_save = zeros(length(cfl_list),length(periodlist),length(amps));

state_iters = -1;
ic_iters = 10;
clearfig = 1;

% semilag_fn = @interp_fourier;
semilag_fn = @interp_cubic;

for aidx = 1:length(amps)
    
    amp_hill = amps(aidx)*phi0;

    for pidx = 1:length(periodlist)
        periods = periodlist(pidx);
        %periods = ceil(Lx/lam_hill); % Number of periods in the full domain (round up)
        k_hill = 2*pi*periods/Lx; % Wavenumber of hill
        hill_fn = @(x) amp_hill*sin((x-dx)*k_hill);

        % Generate initial state and grid
        [phi, u, x_phi, x_u] = gen_state(Lx,Nx,phi0,u0,hill_fn,state_iters);
        hill_phi = hill_fn(x_phi);

        % Get differential operators
        [Dx_phi, Dx_u] = fd_ops(Nx,Lx);
        imat_u_phi = semilag_fn(x_phi,x_u); % Average from u to phi points for semilag
        % Uncomment below for spectral derivatives rather than FD
%         [imat_u_phi, Dx_u] = semilag_fn(x_phi,x_u); % Average from u to phi points for semilag
%         [~, Dx_phi] = semilag_fn(x_u,x_phi);
        


        maxadv = zeros(Ncfl,1);
        maxcn = zeros(Ncfl,1);
        ic_resids = zeros(Ncfl,1);
%         phio = phi; uo = u;
        for idx = 1:Ncfl
            dt = cfl_list(idx)*dx/u0;
            if (ic_iters >= 0)
                % Get IC from O(epsilon) perturbation expansion
                dxhat = dt/dx*sin(0.5*k_hill*dx);
                theta = -k_hill*u0*dt;
                icmat = [-sin(0.5*theta) phi0*dxhat*cos(0.5*theta); 
                         g*dxhat*cos(0.5*theta) -sin(0.5*theta)];
                ichat = icmat \ ...
                         [0;-g*dxhat*cos(0.5*theta)*amp_hill];
                phi_guess = ichat(1)*sin(k_hill*(x_phi-dx))+phi0;
                u_guess = ichat(2)*sin(k_hill*(x_u-dx))+u0;
                
                % We want to detect when this perturbation expansion is not
                % well-defined (i.e. near-singular).  We can do so by
                % checking the determinant of icmat
                
                % For dt->0 and k_hill*dx->0, the matrix becomes
                % [-k*u0*dt/2 phi0*k*dt/2; g*k*dt/2 -k*u0*dt/2]
                % which has determinant (k*dt/2)^2*(u0^2 - c^2)
                
                % That, modified to (u0^2 + c^2), forms the normalization
                % factor -- if the determinant is a small fraction of this
                % then the leading-order IC is not really a small
                % correction to u0/phi0
                if (abs(det(icmat)) < 0.01*k_hill^2*dt^2/4*(c^2 + u0^2))
%                     phio = phi; uo = u; 
                    maxadv(idx) = NaN;
                    maxcn(idx) = NaN;
                    ic_resids(idx) = NaN;
                    continue;
                end
                % Get initial condition
                [phio, uo, ic_resids(idx), imat_u, imat_du, imat_phi, imat_dphi] = gen_ic(phi_guess,u_guess,x_phi,x_u,dt,...
                                    semilag_fn,hill_phi,...
                                    Dx_phi,Dx_u,imat_u_phi,ic_iters,1e-4);
                % If the process did not converge well, try again with
                % regularization
                if (ic_resids(idx) > 1e-2) 
                    [phio, uo, ic_resids(idx), imat_u, imat_du, imat_phi, imat_dphi] = gen_ic(phi_guess,u_guess,x_phi,x_u,dt,...
                        semilag_fn,hill_phi,...
                        Dx_phi,Dx_u,imat_u_phi,ic_iters,1e-4,1);
                end
                                
%                 display(ic_resids(idx))
            end
%             [depart_u,depart_phi] = semilag_traj(uo,imat_u_phi*uo,...
%                                                  x_phi,x_u,dt,@semilag_fn);
%             [imat_u, imat_du] = semilag_fn(depart_u,x_u);
%             [imat_phi, imat_dphi] = semilag_fn(depart_phi,x_phi);
            maxadv(idx) = max(abs(eig(imat_u)));

            cn_op = get_linearized_op(phio,uo,hill_phi,x_phi,x_u,dt,...
                                    imat_phi,imat_dphi,imat_u,imat_du,...
                                    imat_u_phi,Dx_phi,Dx_u);
%             cn_op = get_cn_op(phio,uo,hill_phi,x_phi,x_u,dt,...
%                                     imat_phi,imat_dphi,imat_u,imat_du,...
%                                     imat_u_phi,Dx_phi,Dx_u);
            maxcn(idx) = max(abs(eig(cn_op)));
            if (mod(idx,floor(Ncfl/10)) == 0)
                fprintf('.')
            end
        end
        maxcn_save(:,pidx,aidx) = maxcn;
        maxadv_save(:,pidx,aidx) = maxadv;
        ic_resids_save(:,pidx,aidx) = ic_resids;
        fprintf(' %d %g %g %g\n',periods,max(ic_resids),max(maxcn)-1,max(maxadv)-1);
%         if (length(cfl_list) > 1) 
%             if (clearfig); clf; end
%             semilogy(cfl_list,max(1e-16,maxadv-1),'r-',...
%                      cfl_list,max(1e-16,maxcn-1),'k-');
%     %          semilogy(cfl_list,max(1e-16,maxcn-1),'k-');
%             title(periods)
%             hold on; 
%             axis([cflmin-1e-6 cflmax+1e-6 1e-8 1])
%             drawnow
%         end
        % Save restart file
%         save instab_restart
    end
end

if (length(periodlist) > 1 && Ncfl > 1)
    clf
    rates = maxcn_save.^(1./(cfl_list'*dx/u0));
    %pcolor(periodlist/Nx,cfl_list,rates); shading flat; colorbar;
    pcolor(periodlist,cfl_list,rates); shading flat; colorbar;
    %xlabel('Normalized hill wavenumber (/\Deltax)')
    xlabel('Hill wavenumber');
    ylabel('Base Courant number')
    title('Maximum growth rate')
end
