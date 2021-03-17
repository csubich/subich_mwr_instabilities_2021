% show_spectrum -- show numerical dispersion relation and growth rates for
% a particular [cn_op, imat_u] pairing.  

% The latter describes the waves of the purely advective system (e.g.
% advection of a tracer), which is generally stable except for the case
% where high-frequency waves alias onto one of the same temporal frequency.

% The former case is more varied, since aliasing is possible between
% upstream and downstream-propagating waves.

if (~exist('plotmode','var'))
    plotmode=2;
end

%% Calculate and label eigenvectors
[cn_eigvec, cn_eigval] = eig(cn_op,'vector');
[cn_modes, cn_phases] = label_eigvecs(cn_eigvec,Nx,Lx);
pos_phase = real(cn_phases) > 0;
neg_phase = real(cn_phases) <= 0;

%% Eigenvalue magnitude, CN operator
if (plotmode > 0)
    if (plotmode == 2)
        figure(3)
        subplot(2,1,1)
    elseif (plotmode == 4)
        figure(3)
        subplot(2,2,2)
    end
    plot(1e3*cn_modes(pos_phase),abs(cn_eigval(pos_phase)),'bo','markersize',4);
    hold on
    plot(1e3*cn_modes(neg_phase),abs(cn_eigval(neg_phase)),'r.','markersize',8);
    xlim(1e3*[-pi/dx pi/dx])
    ylim(1+[-1 1]*abs(max(abs(cn_eigval)-1)));
    hold on;
    plot(1e3*[-k_hill k_hill; -k_hill k_hill],...
         min(ylim) + [0.05 0.05;.95 .95]*(max(ylim)-min(ylim)),'k--')
    hold off
    grid on
    title('Eigenvalue magnitude')
    %xlabel('Wavenumber (1/km)')
    ylabel('Magnitude')

    if (plotmode == 2)
        subplot(2,1,2)
    elseif (plotmode == 4)
        subplot(2,2,4)
    end
    % Phase angles, CN operator
    plot(1e3*cn_modes(pos_phase),angle(cn_eigval(pos_phase)),'bo','markersize',4)
    hold on
    plot(1e3*cn_modes(neg_phase),angle(cn_eigval(neg_phase)),'r.','markersize',8);
    axis([1e3*[-pi/dx pi/dx] -3.2 3.2])
    hold on;
end

modespace = linspace(-pi/dx,pi/dx,512); % Continuous range of modes

% Compute analytic (leading order) phase angles, see cn_dispersion.ipynb
c0 = sqrt(phi0*g);
eta = sin(0.5*modespace*dx)*dt/dx;
omega = atan2(2*c0*eta,1-c0^2*eta.^2);

% Compute positive (downstream-propagating waves) phase angles modulo 2pi
pos_angles = mod(pi - omega - u0*dt*modespace,2*pi)-pi;
pos_angles(abs(pos_angles) > 3.1) = nan; % Suppress connection of discontinuities

% Negative phase angles
neg_angles = mod(pi + omega - u0*dt*modespace,2*pi)-pi;
neg_angles(abs(neg_angles) > 3.1) = nan;

if (plotmode > 0)
    plot(1e3*modespace,pos_angles,'b-','linewidth',0.5);
    plot(1e3*modespace,neg_angles,'r-','linewidth',0.5);

    plot(1e3*[-k_hill k_hill; -k_hill k_hill],...
         min(ylim) + [0.05 0.05;.95 .95]*(max(ylim)-min(ylim)),'k--')
    hold off
    grid on

    title('Numerical dispersion relation')
    xlabel('Wavenumber (1/km)')
    ylabel('\omega\Deltat')
end


%% Show most unstable mode
[maxrate, maxidx] = max(abs(cn_eigval));
if (plotmode > 0)
    figure(4);
    ax1 = subplot(2,1,1);
    plot(x_phi/1e3,real(cn_eigvec(1:Nx,maxidx)))
    title('Most unstable \phi perturbation')
    ax2 = subplot(2,1,2);
    plot(x_phi/1e3,hill_phi)
    title('Topography')

    set(ax1,'outerposition',[0 0.35 1 0.65])
    set(ax2,'outerposition',[0 0 1 0.35])
    set(ax1,'xticklabel','')
    grid(ax1,'on')
    xlim(ax1,[-1 1]*Lx/2e3)
    set(ax1,'yticklabel','')
    set(ax2,'yticklabel','')
    xlabel(ax2,'x (km)')
    xlim(ax2,[-1 1]*Lx/2e3)
end


%% Advection operator
[adv_eigvec, adv_eigval] = eig(imat_u,'vector');
[adv_modes] = label_eigvecs(adv_eigvec,Nx,Lx);

if (plotmode > 0)
    if (plotmode == 2)
        figure(2)
        subplot(2,1,1)
    elseif (plotmode == 4)
        figure(3)
        subplot(2,2,1)
    end
    % Magnitude, advective operator
    plot(1e3*adv_modes,abs(adv_eigval),'k.','markersize',8);
    xlim(1e3*[-pi/dx pi/dx])
    hold on;
    plot(1e3*[-k_hill k_hill; -k_hill k_hill],...
         min(ylim) + [0.05 0.05;.95 .95]*(max(ylim)-min(ylim)),'k--')
    hold off
    grid on
    title('Eigenvalue magnitude (advection)')
    %xlabel('Wavenumber (1/km)')
    ylabel('Magnitude')

    if (plotmode == 2)
        subplot(2,1,2)
    elseif (plotmode == 4)
        subplot(2,2,3)
    end
    % Phase angles, advective operator
    plot(1e3*adv_modes,angle(adv_eigval),'k.','markersize',8);
    axis([1e3*[-pi/dx pi/dx] -3.2 3.2])

    hold on;
end

adv_angles = mod(pi - u0*dt*modespace,2*pi)-pi;
adv_angles(abs(adv_angles) > 3.1) = nan;

if (plotmode > 0)
    plot(1e3*modespace,adv_angles,'k-')


    plot(1e3*[-k_hill k_hill; -k_hill k_hill],...
         min(ylim) + [0.05 0.05;.95 .95]*(max(ylim)-min(ylim)),'k--')
    hold off
    grid on
    title('Numerical dispersion relation')
    xlabel('Wavenumber (1/km)')
    ylabel('\omega\Deltat')
end