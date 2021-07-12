% Show interesting spectra of various instability cases

% Basic paramters -- copied from paper_instab.m


dx = 1000;
if (~exist('Nx','var'))
    Nx = 513;
end
%Nx = 129;
%Nx = 513;
Lx = dx*Nx;
u0 = 20; % Mean velocity
phi0 = 6000; % Mean water depth
global g
g = 9.81;
c = sqrt(g*phi0); % phase speed

% Grid
% phi points are located one half-cell away from the periodic boundary
x_phi = Lx*(-0.5 + (-0.5 + (1:Nx)')/Nx);
x_phi = x_phi(:);
% u-points are located one half-cell to the right of phi-points, s.t.
% x_u(1) - x_phi(1) = dx/2.  x_u(end) is the periodic boundary at Lx/2.
x_u = x_phi + dx/2;
x_u = x_u(:);

% Finite difference operators
[Dx_phi, Dx_u] = fd_ops(Nx,Lx);
imat_u_phi = interp_fourier(x_phi,x_u); % Average from u to phi points for semilag

% Hill amplitude
% 1.67e-2 -> 100m hills on 6km phi0
if (~exist('amp_mult','var'))
    amp_mult = 1;
end
amp_hill = amp_mult*1.667e-2*phi0;

if (~exist('CASE','var'))
    CASE = 1;
end

%% CASE 1

if CASE == 1
% This case features a purely advective instability.  The semi-Lagrangian
% operator itself has eigenvalues of magnitude > 1, where high-frequency
% waves alias across the Nyquist limit onto an aliased mode with the
% same temporal frequency.  This is only possible when the cfl > 1.

%wavelength = 257e3/39;
wavelength = 10*dx;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
%cfl = 1.18;
cfl = 1.11;
ic_iters = 12;
end

%% CASE 2 -- Self-resonance
if CASE == 2
% This case features an orography-induced resonance.  The group velocity
% corresponding to the upstream-propagating wave reaches 0 and reverses
% sign.  As a result, the upstream-propagating wave around k=+-1.1617e-6
% shares a temporal frequency with the downstream-propagating wave around
% k=+-(k_crit + k_hill).  The two separate waves combine into a
% growing/decaying pair, leading to an unstable mode.

%wavelength = 257e3/20;
%wavelength = 513e3/40;
wavelength = 513e3/117;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
%cfl = 0.9224;
%cfl = 0.8638;
cfl = 0.6653;

% This also happens at higher Courant numbers, see parameters:
% wavelength = 257e3/19
% wavelength = Lx / ceil(Lx/wavelength);
% k_hill = 2*pi/wavelength; 
% cfl = 2.038;
% dt = cfl*dx/u0;

ic_iters = 12;
end

%% CASE 3 -- Wave Nyquist Resonance
if CASE == 3
% This case is somewhat analaogous to case 1.  Here, a wave of spatial
% frequency k_crit gains a side-lobe at frequency k_crit + k_hill, which
% aliases onto a wave (k_crit + k_hill - 2*k_nyq) that has the same
% frequency.  This is visible here for the downstream-propagating wave near
% k=+-2e-3; the aliased portion of the wave is visible as the missing
% wavenumbers in the dispersion relation.

% wavelength = 257e3/57
% wavelength = Lx / ceil(Lx/wavelength);
% k_hill = 2*pi/wavelength; 
% cfl = 0.266;
% dt = cfl*dx/u0;

% Alternative parameters:
wavelength = 513e3/76;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
cfl = 0.2485;
ic_iters = 12;

end

%% CASE 4 -- Temporal-alising resonance
% This case features resonance between the downstream-propagating wave at
% low wavenumber and its higher-wavenumber counterpart.

% Here, the positive-going wave at about k*=+-7.33e-5 (wavelength 8.57km)
% has nearly the same frequency as that at +-7.33e-5 -+ k_hill (wavelength
% 5.47km) (modulo 2pi).  This aliasing is a feature of the high CFL number
% (1.683), which combines with the relatively fast gravity-wave speed to
% to give a dispersion relation that wraps around at a low weavenumber.

% This instability does not directly result from spatial aliasing -- all
% the relevant waves are well-resolved on the grid.  Additionally, it does
% not require spurious orographic resonance, since the group velocity of
% this wave is bounded well away from 0 for the full range of spatial
% frequencies.  (This is not the case for the negative-going wave, and an
% instability of case 2 is evident at lower amplitude).
if CASE == 4
wavelength = 257e3/50;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
% cfl = 1.683;
cfl = 1.672;
ic_iters = 12;
% ic_iters = 0;
end
%% CASE 5 -- Strong spurious resonance
% This case is much like case 4, only now it's the upstream-propagating
% branch of waves that has aliasing.  Here, the unstable mode is
% k=+-2.45e-5 (mode 1, Lx=257km), which aliases onto very-nearly the
% topographic frequency.  

% The negative-going wave branch does not cross omega=+-pi, so this 
% aliasing arises from the numerically-induced dispersion
if CASE == 5
wavelength = 257e3/50;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
cfl = 3.5;
ic_iters = 12;
end

%% Working case (parameter-space exploration)
if 0
wavelength = 513e3/123;
wavelength = Lx / ceil(Lx/wavelength);
k_hill = 2*pi/wavelength; 
cfl = cfl_list(157) %1.737;
ic_iters = 5;
end


dt = cfl*dx/u0;

hill_fn = @(x) amp_hill*sin((x-dx)*k_hill);
hill_phi = hill_fn(x_phi);
% Generate initial condition
dxhat = dt/dx*sin(0.5*k_hill*dx);
theta = -k_hill*u0*dt;
icmat = [-sin(0.5*theta) phi0*dxhat*cos(0.5*theta); 
         g*dxhat*cos(0.5*theta) -sin(0.5*theta)];
ichat = icmat \ ...
         [0;-g*dxhat*cos(0.5*theta)*amp_hill];
phi_guess = ichat(1)*sin(k_hill*(x_phi-dx))+phi0;
u_guess = ichat(2)*sin(k_hill*(x_u-dx))+u0;

% Revise initial condition
[phio, uo, ic_resid, imat_u, imat_du, imat_phi, imat_dphi] = gen_ic(phi_guess,u_guess,x_phi,x_u,dt,...
                    @interp_fourier,hill_phi,...
                    Dx_phi,Dx_u,imat_u_phi,ic_iters,1e-4);
if (ic_resid > 1e-2 && ic_iters > 0) 
    [phio, uo, ic_resids(idx), imat_u, imat_du, imat_phi, imat_dphi] = gen_ic(phi_guess,u_guess,x_phi,x_u,dt,...
        @interp_fourier,hill_phi,...
        Dx_phi,Dx_u,imat_u_phi,ic_iters,1e-4,1);
end
display(ic_resid)

% Get operator
cn_op = get_linearized_op(phio,uo,hill_phi,x_phi,x_u,dt,...
                        imat_phi,imat_dphi,imat_u,imat_du,...
                        imat_u_phi,Dx_phi,Dx_u);
                    
show_spectrum

% Display interesting info
fprintf('Hill wavenumber %.3g, %.2f dx\n',k_hill,2*pi/(k_hill*dx));
fprintf('Advective Courant number %.2f\n',dt*u0/dx);
[maxrate, maxidx] = max(abs(cn_eigval));
fprintf('Maximum wave growth rate: %.3f / step\n',maxrate);
fprintf('Unstable mode is approximately %.3e\n',cn_modes(maxidx))
% 
% figure(4);
% ax1 = subplot(2,1,1);
% plot(x_phi/1e3,real(cn_eigvec(1:Nx,maxidx)))
% title('Most unstable \phi perturbation')
% ax2 = subplot(2,1,2);
% plot(x_phi/1e3,hill_phi)
% title('Topography')
% 
% set(ax1,'outerposition',[0 0.35 1 0.65])
% set(ax2,'outerposition',[0 0 1 0.35])
% set(ax1,'xticklabel','')
% grid(ax1,'on')
% xlim(ax1,[-1 1]*Lx/2e3)
% set(ax1,'yticklabel','')
% set(ax2,'yticklabel','')
% xlabel(ax2,'x (km)')
% xlim(ax2,[-1 1]*Lx/2e3)