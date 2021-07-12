% Script to demonstrate general accuracy of the time-stepper, with a stable
% versus unstable case (against controls)

% Domain and problem parameters
dx = 1000;
Nx = 513;
Lx = dx*Nx;
u0 = 20;
phi0 = 6000;
amp_hill = 100;
global g
g = 9.81;
c = sqrt(g*phi0);

% Staggered Grid
x_phi = Lx*(-0.5 + (-0.5 + (1:Nx)')/Nx);
x_phi = x_phi(:);
% u-points are located one half-cell to the right of phi-points, s.t.
% x_u(1) - x_phi(1) = dx/2.  x_u(end) is the periodic boundary at Lx/2.
x_u = x_phi + dx/2;
x_u = x_u(:);

% Hill -- case 4 of interesting_spectra
hill_wavelength = 5*dx;
hill_wavelength = Lx/ ceil(Lx/hill_wavelength);
k_hill = 2*pi/hill_wavelength;

hill_fn = @(x) amp_hill*sin((x-dx)*k_hill);
hill_phi = hill_fn(x_phi);
hill_u = hill_fn(x_u);

% First-order steady state flow
phio = phi0*(1 - g*hill_phi./(phi0*g - u0.^2));
uo = u0*(1 + g*hill_u./(phi0*g - u0.^2));

% add 1m perturbation to phi
rng(20200707);
phio = phio + randn(Nx,1);
pert_eps = 0; % Clear builtin perturbation in timestep

% Case 1: high CFL
cfl = 1.5;
dt = cfl*dx/u0;
clear dtlist;
timestep;
phia_1p5 = phia;

% Case 2: low CFL
cfl = 0.5;
dt = cfl*dx/u0;
clear dtlist;
timestep;
phia_0p5 = phia;

% Reference case, with ode45 and spectral method
kk = [0:floor(Nx/2) -floor(Nx/2):-1]'*2*pi/Lx;
uo_phi = u0*(1 + g*hill_phi./(phi0+g - u0.^2)); % Get u on phi grid

fftd = @(f) ifft(1i.*kk.*fft(f));
sw_rhs = @(t,f) [-fftd(f(1:Nx).*f(Nx+(1:Nx)));-fftd(0.5*f(Nx+(1:Nx)).^2 + g*(hill_phi + f(1:Nx)))];

% Set tiny tolerances for ode45
ode_opts  = odeset('AbsTol',1e-7,'RelTol',1e-7);
% Evaluate
[tt, yy] = ode45(sw_rhs,[0 3600*4],[phio;uo],ode_opts);

% Plot

subplot(3,1,1); plot(x_phi/1e3,yy(end,1:Nx)'+hill_phi);
title('ode45 reference')
ylabel('\phi + H')
set(gca,'xticklabels',[]);
xlim([-Lx/2,Lx/2]/1e3);
subplot(3,1,2); plot(x_phi/1e3,phia_0p5+hill_phi);
title('Semi-Lagrangian, C=0.5');
ylabel('\phi + H')
set(gca,'xticklabels',[]);
xlim([-Lx/2,Lx/2]/1e3);
subplot(3,1,3); plot(x_phi/1e3,phia_1p5+hill_phi);
title('Semi-Lagrangian, C=1.5');
ylabel('\phi + H')
xlim([-Lx/2,Lx/2]/1e3);
xlabel('x (km)');

fprintf('ODE45 mean timestep %f, C=%.2e\n',mean(diff(tt)),mean(diff(tt))*u0/dx);