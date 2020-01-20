function plot3_pretty(Parameters,Options,t,y,m,n,sij,spikes,raster,rr,Ihistory,sijhistory)

N = Parameters.n1 + Parameters.n2;
tauavg = Parameters.tauavg;
dt = Parameters.dt;
tmax = Parameters.tmax;

nplots = 8;

tstart = dt;
tend = tmax;

tw = tstart/dt:tend/dt;

fig = figure('Position',[800 500 1500 1200]);

ax1 = subplot(nplots,1,1);
hold on
plot(t(tw)/1000,(spikes(tw)./tauavg)*1000)
title('Network Activity (Spikes/Neuron/Second)')
xlim([5 35])

%%%%%%%%%%%%%%%%% plot 2 - raster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = subplot(nplots,1,2); hold on;

for i = 1:N
    plot(t/1000,raster(:,i)*i,'.k'); hold on;
end
    
ylim([1 N])
xlim([5 35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title('Raster Plot')


ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,y(tw,rr))
title(['Example neuron`s synaptic depression (y)'])
xlim([5 35])

ax5 = subplot(nplots,1,4);
plot(t(tw)/1000,m(tw,rr))
title(['Example random neuron`s fast synaptic depression term (m)'])
xlim([5 35])

ax6 = subplot(nplots,1,5);
plot(t(tw)/1000,n(tw,rr))
title(['Example neuron`s slow synaptic depression term (n)'])
xlim([5 35])

ax4 = subplot(nplots,1,6);
plot(t(tw)/1000,sij(tw,rr))
title(['Example neuron`s synaptic current output (s)'])
xlim([5 35])

ax7 = subplot(nplots,1,7);
plot(t(tw)/1000,sijhistory(tw))
title(['Example neuron`s synaptic current input'])
xlim([5 35])


ax8 = subplot(nplots,1,8);
plot(t(tw)/1000,Ihistory(tw))
xlabel('Time (s)')
title(['Example neuron`s total applied current (I+delta*syncurrent*syndepression)'])
xlim([5 35])

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'x')

saveas(gcf,'zoomin.png')

%saveas(fig,['dataset' mat2str(length(dir('*.png'))+1) '.png'])
