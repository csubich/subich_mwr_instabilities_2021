% show_spectrum -- show numerical dispersion relation and growth rates for
% a particular [cn_op, imat_u] pairing.  

% The latter describes the waves of the purely advective system (e.g.
% advection of a tracer), which is generally stable except for the case
% where high-frequency waves alias onto one of the same temporal frequency.

% The former case is more varied, since aliasing is possible between
% upstream and downstream-propagating waves.

plot_stride = 10;
percentile = 0.50;

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
    % Plot most unstable points
    med_eigval = median(abs(cn_eigval));
    range = max(max(abs(cn_eigval)) - med_eigval, med_eigval - min(abs(cn_eigval)));
    max_val = med_eigval + percentile*range;
    min_val = med_eigval - percentile*range;
    
    big_pos = pos_phase' & ((abs(cn_eigval) > max_val) | (abs(cn_eigval) < min_val));
    big_neg = neg_phase' & ((abs(cn_eigval) > max_val) | (abs(cn_eigval) < min_val));
    
    % Plot the 10% extreme values and 1/5 of the remainder
    plot(1e3*cn_modes(big_pos),abs(cn_eigval(big_pos)),'bo','markersize',4);
    
    hold on
    plot(1e3*cn_modes(big_neg),abs(cn_eigval(big_neg)),'r.','markersize',8)
    
    % Sort modes to plot the remainder
    pos_idx = find(pos_phase);
    [~,ii] = sort(cn_modes(pos_idx));
    pos_idx = pos_idx(ii);
    
    neg_idx = find(neg_phase);
    [~,ii] = sort(cn_modes(neg_idx));
    neg_idx = neg_idx(ii);
    
    plot(1e3*cn_modes(pos_idx(1:plot_stride:end)),...
        abs(cn_eigval(pos_idx(1:plot_stride:end))),'bo','markersize',4);
    plot(1e3*cn_modes(neg_idx(1:plot_stride:end)),...
        abs(cn_eigval(neg_idx(1:plot_stride:end))),'r.','markersize',8);
    xlim(1e3*[-pi/dx pi/dx])
    ylim(1+[-1 1]*abs(max(abs(cn_eigval)-1)));
    
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
    plot(1e3*cn_modes(big_pos),angle(cn_eigval(big_pos)),'bo','markersize',4);
    hold on
    plot(1e3*cn_modes(pos_idx(1:plot_stride:end)),...
        angle(cn_eigval(pos_idx(1:plot_stride:end))),'bo','markersize',4)
    plot(1e3*cn_modes(big_neg),angle(cn_eigval(big_neg)),'r.','markersize',8);
    plot(1e3*cn_modes(neg_idx(1:plot_stride:end)),...
        angle(cn_eigval(neg_idx(1:plot_stride:end))),'r.','markersize',8);
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
    modenorm = max(abs(cn_eigvec(1:Nx,maxidx)))
    plot(x_phi/1e3,real(cn_eigvec(1:Nx,maxidx))/modenorm)
    title('Most unstable \phi perturbation')
    ax2 = subplot(2,1,2);
    plot(x_phi/1e3,hill_phi)
    title('Topography')

    set(ax1,'outerposition',[0 0.35 1 0.65])
    set(ax2,'outerposition',[0 0 1 0.35])
    set(ax1,'xticklabel','')
    grid(ax1,'on')
    xlim(ax1,[-1 1]*Lx/2e3)
    ylabel(ax1,'\phi'' (m)')
    %set(ax1,'yticklabel','')
    %set(ax2,'yticklabel','')
    xlabel(ax2,'x (km)')
    ylabel(ax2,'H (m)')
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

    med_eigval = median(abs(adv_eigval));
    range = max(max(abs(adv_eigval)) - med_eigval, med_eigval - min(abs(adv_eigval)));
    max_val = med_eigval + percentile*range;
    min_val = med_eigval - percentile*range;
    
    big_adv = ((abs(adv_eigval) > max_val) | (abs(adv_eigval) < min_val));
    
    
    [~,adv_idx] = sort(adv_modes);
    
    
    plot(1e3*adv_modes(big_adv),abs(adv_eigval(big_adv)),'k.','markersize',8)
    hold on;
    plot(1e3*adv_modes(adv_idx(1:plot_stride:end)),...
             abs(adv_eigval(adv_idx(1:plot_stride:end))),'k.','markersize',8);
    xlim(1e3*[-pi/dx pi/dx])
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
    plot(1e3*adv_modes(big_adv),angle(adv_eigval(big_adv)),'k.','markersize',8);
    hold on
    plot(1e3*adv_modes(adv_idx(1:plot_stride:end)),...
             angle(adv_eigval(adv_idx(1:plot_stride:end))),'k.','markersize',8);
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