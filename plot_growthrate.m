% Plot growth rate from run of paper_instab

clf
[foo,bar] = contourf(periodlist*dx/Lx,cfl_list,3600*log(rates)/log(10),linspace(0.17,4,16));
hold on
[qux quux] = contourf(periodlist*dx/Lx,cfl_list,isnan(rates),[0.1 1]);
hold off

set(bar,'edgecolor','none')
set(quux,'edgecolor',0.75*[1 1 1],'facecolor',0.75*[1 1 1])
cb = colorbar;
caxis([0 2])
colormap(parula(16))
set(cb,'ylim',[0.17 2]);
set(cb,'ytick',[0.1761, 0.477 1 1.477 2])
set(cb,'ticklabels',{'1.5','3','10','30','100'})
%ylabel(cb,'Growth factor per hour')

fprintf('Global maximum growth rate %.2e\n',10^max(3600*log(rates(:))/log(10)))