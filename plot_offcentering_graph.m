% Present waterfall plot of growth rates versus off-centering parameter and
% hill amplitude, using existing test_with_offcentering script and
% associated magic parameters.  

% Case paramters are:
% cfl = 1.683;
% hill_width = 257e3/50;

% Multipliers on the base hill amplitude: 0.1 -> 10 (10m -> 1km amplitude)

if (~exist('plot_paper_figs','var')) % Plotting figures for paper or talk?
    plot_paper_figs = 1;
end
if (~exist('plot_talk_figs','var'))
    plot_talk_figs = 1;
end
%hill_mults = [0.1;0.225;0.375;0.5;0.75;1;2.25;3.75]%;5;7.5;10];
%hill_amps = [10 25 37.5 50 75 100 250 375];
%hill_amps = [10 20 40 60 80 100 200 400];
hill_amps = logspace(1,3,11);
mark1 = 1; mark2 = 6; mark3 = 9;
% hill_mults = hill_amps/phi0/1.67e-2;
num_mults = length(hill_amps);

% Off-centering paramters
alphas = linspace(0.5,0.6,81);
num_alphas = length(alphas);

maxeig = zeros(num_alphas,num_mults);

for idx1 = 1:num_alphas
    alpha = alphas(idx1);
    for idx2 = 1:num_mults
        %amp_mult = hill_mults(idx2);
        amp_hill = hill_amps(idx2);
        test_with_offcentering;
        maxeig(idx1,idx2) = max(abs(eig(cn_me)));
        fprintf('%.3f %7.2f %.2e\n',alpha,amp_hill,maxeig(idx1,idx2)-1);
    end
end

amp_factor = (maxeig.^(3600/dt)-1);
        
clf;
plot(alphas,amp_factor,'k-');
set(gca,'yscale','log');
%axis([0.5 0.6 1e-6 2])
axis([0.5 0.6 1e-5 1.1e4])
ylabel('Amplification (/hr)');
xlabel('Off-centering \alpha');

plot_alpha = floor(num_alphas/2);
text(alphas(plot_alpha),0.95*amp_factor(plot_alpha,mark1),...
    sprintf('A_{hill}=%.0fm',hill_amps(mark1)),...
    'verticalalignment','top')

text(alphas(plot_alpha)-0.027,0.98*amp_factor(plot_alpha,mark2),...
    sprintf('A_{hill}=%.0fm',hill_amps(mark2)),...
    'verticalalignment','bottom')

text(alphas(plot_alpha),1.05*amp_factor(plot_alpha,mark3),...
    sprintf('A_{hill}=%.0fm',hill_amps(mark3)),...
    'verticalalignment','bottom')

if (plot_talk_figs)
%     title(sprintf('Growth rates, %.0f\\Deltax & C=%.2f',(2*pi/k_hill)/dx,dt*u0/dx));
%     set(gcf,'papersize',0.75*[5,3.8]);
%     set(gcf,'paperposition',[0 0 5 3.8]*0.75);
% 
%     print -dpdf -painters talk_figs/fig_alphas.pdf

    set(gca,'ytick',[1e-4,1e-2,1,1e2,1e4])
    set(gca,'yticklabels',{'1+10^{-4}','1+10^{-2}','2','10^{2}','10^{4}'})
    title('Off-centering growth rates')
    set(gcf,'papersize',0.75*[5,3.8]);
    set(gcf,'paperposition',[0 0 5 3.8]*0.75);
    print -dpdf -painters pdes_figs/fig_offc_alphas.pdf
end
if (plot_paper_figs)
    set(gca,'ytick',[1e-4,1e-2,1,1e2,1e4])
    set(gca,'yticklabels',{'1+10^{-4}','1+10^{-2}','2','10^{2}','10^{4}'})
    title('Off-centering growth rates')
    set(gcf,'papersize',[5 3],'paperposition',[0 0 5 3]);
    print -dpdf -painters paper_figs/fig_offc_alphas.pdf
end
    