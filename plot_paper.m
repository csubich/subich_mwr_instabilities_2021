% Consolidated script to generate figures for the to-be-submitted paper

% Default figure size/position settings
%close all
% 5in by 3in figures
set(0,'defaultfigurepapersize',[5 3]);
clear

GROWTHRATE_PLOTS = 0;
TIME_PLOTS = 1;
ALT_SEMILAG_PLOTS = 0;
OFFC_PLOTS = 0;
TIME_CONVERGENCE = 0;

%% Growth rate figure
if GROWTHRATE_PLOTS

% This will be assembled in parts to update labels for various unstable
% regions "live", consistently with subsequent figures

load growthrate_fourier_512
figure(1);
plot_growthrate
% Set on-screen position of this figure because the arrows are plotted in
% normalized coordinates; a change of aspect ratio at print time upsets the
% figure.
set(gcf,'paperposition',[0 0 5 3])
set(gcf,'position',[0 0 500 300]);

set(gca,'xtick',[0.05 0.1 0.2 0.25 0.3333])
set(gca,'xticklabels',{'20\Deltax','10\Deltax','5\Deltax','4\Deltax','3\Deltax'})
xlabel('Topography wavelength')
ylabel('Courant number')
title('Maximum amplification per hour');


% Draw arrow, from https://www.mathworks.com/matlabcentral/answers/346297-how-to-draw-an-arrow-using-non-normalized-coordinates
% Difficult because arrows are drawn in normalized figure units, not data
% units
pos = get(gca,'position');
arrow_x = [0.27 0.25];
arrow_y = [1.9 2];
hold on
plot(arrow_x(2),arrow_y(2),'k.','markersize',28);
plot(arrow_x(2),arrow_y(2),'r.','markersize',20);
rx=xlim(gca);ry=ylim(gca);
cx=pos(3)/diff(rx);cy=pos(4)/diff(ry)
annotation('textarrow',pos(1)+cx*(arrow_x-rx(1)),pos(2)+cy*(arrow_y-ry(1)),...
    'String','(a)','textmargin',0.01)

%% Advective spatial/temporal aliasing
%Nx = 129; % Don't need many points
Nx = 513;
plotmode = 4; % Advective instability in a 2x2 layout
CASE=1; interesting_spectra

figure(3);
set(gcf,'paperposition',[0 0 5 3])
print -dpdf -painters paper_figs/fig_advect_spectrum.pdf
% Add marker to growth rate graph
figure(1);
hold on;
arrow_x = k_hill*dx/(2*pi) + [-0.02,0];
arrow_y = cfl + [-0.1;0];
% text(k_hill*dx/(2*pi),cfl,'(b)','verticalalign','top')
plot(arrow_x(2),arrow_y(2),'k.','markersize',28);
plot(arrow_x(2),arrow_y(2),'r.','markersize',20);
rx=xlim(gca);ry=ylim(gca);
cx=pos(3)/diff(rx);cy=pos(4)/diff(ry)
annotation('textarrow',pos(1)+cx*(arrow_x-rx(1)),pos(2)+cy*(arrow_y-ry(1)),...
    'String','(b)','textmargin',0.01)

%% Wave spatial aliasing
%Nx = 127;
Nx = 513;
CASE=3; 
plotmode=2;  % 1x2 layout -- no need to plot advection
interesting_spectra

figure(3)
set(gcf,'paperposition',[0 0 5 3])
print -dpdf -painters paper_figs/fig_wave_space_spectrum.pdf

figure(1);
hold on;
% text(k_hill*dx/(2*pi),cfl,'(c)','verticalalign','bottom')

arrow_x = k_hill*dx/(2*pi) + [-0.02,0];
arrow_y = cfl + [.2;0];
% text(k_hill*dx/(2*pi),cfl,'(b)','verticalalign','top')
plot(arrow_x(2),arrow_y(2),'k.','markersize',28);
plot(arrow_x(2),arrow_y(2),'r.','markersize',20);
rx=xlim(gca);ry=ylim(gca);
cx=pos(3)/diff(rx);cy=pos(4)/diff(ry)
annotation('textarrow',pos(1)+cx*(arrow_x-rx(1)),pos(2)+cy*(arrow_y-ry(1)),...
    'String','(c)','textmargin',0.01)

%% Wave temporal aliasing
% Two cases
%Nx = 257;
Nx = 513;
CASE=4;
plotmode=2;
interesting_spectra

figure(3)
set(gcf,'paperposition',[0 0 5 3])
print -dpdf -painters paper_figs/fig_wave_timea_spectrum.pdf
figure(4)
set(gcf,'paperposition',[0 0 5 3])
print -dpdf -painters paper_figs/fig_wave_timea_profile.pdf

figure(1)
hold on;
% text(k_hill*dx/(2*pi),cfl,'(d)','verticalalign','top')
arrow_x = k_hill*dx/(2*pi) + [0.02,0];
arrow_y = cfl + [-0.1;0];
% text(k_hill*dx/(2*pi),cfl,'(b)','verticalalign','top')
plot(arrow_x(2),arrow_y(2),'k.','markersize',28);
plot(arrow_x(2),arrow_y(2),'r.','markersize',20);
rx=xlim(gca);ry=ylim(gca);
cx=pos(3)/diff(rx);cy=pos(4)/diff(ry)
annotation('textarrow',pos(1)+cx*(arrow_x-rx(1)),pos(2)+cy*(arrow_y-ry(1)),...
    'String','(d)','textmargin',0.01)

% Finally, print growth-rate figure
figure(1)
print -dpdf -painters paper_figs/fig_growthrate.pdf

end

%% Off-centering
if (OFFC_PLOTS)
%Nx = 257;
Nx = 513;
CASE=4;
plotmode=0;
interesting_spectra; close all

% Plot single spectrum with off-centering
alpha = 0.525;
test_with_offcentering; cn_op = cn_me; plotmode = 2; show_spectrum;
% Spectrum plot is in figure 3, update axis limits to show the eigenvalues
figure(3); subplot(2,1,1);
ylim([0.9 1.0028])
hold on;
plot(1e3*[-k_hill k_hill; -k_hill k_hill],...
     min(ylim) + [0.05 0.05;.95 .95]*(max(ylim)-min(ylim)),'k--')
hold off
fprintf('Maximum eigenvalue (alpha = %f) is 1+%.6e\n',alpha,max(abs(cn_eigval))-1)

set(gcf,'paperposition',[0 0 5 3])
print -dpdf -painters paper_figs/fig_offc_spectrum.pdf


plot_paper_figs = 1; plot_talk_figs = 0;
plot_offcentering_graph; % Includes own print statement

% % With off-centering
% 
% alpha = 0.52; amp_mult = 1.0; test_with_offcentering
% cn_op = cn_me; plotmode = 0; show_spectrum;
% figure(3)
% %print -dpdf -painters paper_figs/fig_a52_m1_spectrum.pdf
% 
% 
% alpha = 0.55; amp_mult = 4.0; test_with_offcentering
% cn_op = cn_me; plotmode = 0; show_spectrum;
% figure(3)
% %print -dpdf -painters paper_figs/fig_a55_m4_spectrum.pdf

end

%% Timestepping
if TIME_PLOTS
plotmode = 0; % No plotting spectra
Nx = 513;
CASE=4; interesting_spectra; close all

% Baseline case: random perturbation, full nonlinear solver, no timestep
% shenanigans
fprintf('Timestepping: basic case\n');

INIT='random'; SOLVER='newton';
FINTIME=8;
iter_nonlin_limit = 10;
% Disable diagnostic prints/plots
print_stats = 0;
ploteach=nan;

    
timestep

% Plot results
%close all
figure(1);clf

subplot(2,1,1)
plot(x_phi/1e3,phia + hill_phi,'-'); hold on; plot(x_phi/1e3,phio + hill_phi,'--')
xlabel('x (km)')
ylabel('\phi + H (m)')
xlim(1e-3*[-0.5 0.5]*Lx);
title(sprintf('Layer height, t=%.0fh',times(end)/3600))

subplot(2,1,2)
semilogy(times/3600,phi_amp,'linewidth',1.5)
dasht = 3600*linspace(2,4,32);
dashamp = 1e0*max(abs(cn_eigval)).^((dasht - min(dasht))/dt);
hold on;
semilogy(dasht/3600,dashamp,'k--','linewidth',0.5)
hold on
xlabel('t (h)')
ylabel('Amplitude (m)')
title('Perturbation RMS amplitude')
xlim([0,FINTIME]);
text(dasht(16)/3600,dashamp(16),'Theoretical rate','horizontalalignment','left','verticalalignment','top')
set(gca,'ytick',[1 10 100],'yticklabels',{'1','10','100'})

set(gcf,'paperposition',[0 0 5 3])

drawnow

print -dpdf -painters paper_figs/fig_timestep_base.pdf

times_0 = times;
amps_0 = phi_amp;


% Timestepping: timestep scheduling
% fprintf('Timestep scheduling\n')
% INIT='random'; SOLVER='newton';
% iter_nonlin_limit = 10;
% FINTIME=8;
% 
% dtlist = [0.5 0.5 0.5 0.5 3]*dt;
% timestep;
% 
% times_hop = times;
% amps_hop = phi_amp;
% 
% figure(2); clf
% subplot(2,1,1)
% semilogy(times_hop/3600,amps_hop,'-','linewidth',1.5);
% hold on
% semilogy(times_0/3600,amps_0,'k--','linewidth',0.5);
% title('Perturbation growth - timestep variation')
% ylabel('Amplitude')
% 
% 
% xlim([0,FINTIME]);
% set(gca,'ytick',[1 10 100],'yticklabels',{'1','10','100'})
% set(gca,'xticklabels',[])


% Semi-implicit timestepping

fprintf('Truncated iterations\n')
SOLVER='linear'; iter_nonlin_limit = 2;
FINTIME=8; dtlist = dt;

timestep
times_ilin = times;
amps_ilin = phi_amp;

figure(2)
subplot(2,1,1)
semilogy(times_ilin/3600,amps_ilin,'-','linewidth',1.5);
hold on
semilogy(times_0/3600,amps_0,'k--','linewidth',0.5);
title('Perturbation growth - iterative solve (2 iterations)')
ylabel('Amplitude (m)')
xlabel('Time (h)')


xlim([0,FINTIME]);
ylim([0.1, max(amps_0)]);
set(gca,'ytick',[0.1 1 10 100],'yticklabels',{'0.1','1','10','100'})

fprintf('Iterative solve (1 iteration / semi-implicit)\n')
SOLVER='linear'; iter_nonlin_limit=1;
FINTIME=8; dtlist = dt;

timestep
times_slin = times;
amps_slin = phi_amp;

figure(2)

subplot(2,1,2)
semilogy(times_slin/3600,amps_slin,'-','linewidth',1.5);
hold on
semilogy(times_0/3600,amps_0,'k--','linewidth',0.5);
title('Iterative solve (1 iteration / semi-implicit)')
ylabel('Amplitude (m)')
xlabel('Time (h)')


xlim([0,FINTIME]);
ylim([0.1, max(amps_0)]);
set(gca,'ytick',[0.1 1 10 100],'yticklabels',{'0.1','1','10','100'})


set(gcf,'paperposition',[0 0 5 3])
print -painters -dpdf paper_figs/fig_timestep_alt.pdf
end

%% Alternative semilag operators
if ALT_SEMILAG_PLOTS
    % Linear
    load growthrate_linear_512
    figure(1);
    plot_growthrate

    set(gca,'xtick',[0.05 0.1 0.2 0.25 0.3333])
    set(gca,'xticklabels',{'20\Deltax','10\Deltax','5\Deltax','4\Deltax','3\Deltax'})
    xlabel('Topography wavelength')
    ylabel('Courant number')
    title('Amplification factor (linear)');
    set(gcf,'paperposition',[0 0 5 3])
    print -dpdf -painters paper_figs/fig_growth_linear.pdf
    
    load growthrate_cubic_512
    figure(1);
        plot_growthrate

    set(gca,'xtick',[0.05 0.1 0.2 0.25 0.3333])
    set(gca,'xticklabels',{'20\Deltax','10\Deltax','5\Deltax','4\Deltax','3\Deltax'})
    xlabel('Topography wavelength')
    ylabel('Courant number')
    title('Amplification factor (cubic)');
    set(gcf,'paperposition',[0 0 5 3])
    print -dpdf -painters paper_figs/fig_growth_cubic.pdf
        
end

if TIME_CONVERGENCE
    clear;
    time_convergence;
    set(gcf,'paperposition',[0 0 5 4]);
    set(gcf,'papersize',[5 4]);
    set(gcf,'position',[0 0 500 400])
    print -dpdf -painters paper_figs/fig_time_convergence.pdf
end
    
    