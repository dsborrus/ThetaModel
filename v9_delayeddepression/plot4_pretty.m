function plot4_pretty(Parameters,Options,t,y,m,n,sij,spikes,raster,rr,Ihistory,sijhistory)

N = Parameters.n1 + Parameters.n2;
tauavg = Parameters.tauavg;
dt = Parameters.dt;
tmax = Parameters.tmax;

nplots = 4;

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
plot(t(tw)/1000,mean(y(tw,:),2))
title(['Mean synaptic depression'])
xlim([5 35])

ax4 = subplot(nplots,1,4);
plot(t(tw)/1000,mean(sij(tw,:),2))
title(['Mean synaptic current output (s)'])
xlim([5 35])
xlabel('Time (s)')

saveas(gcf,'zoomout.png')

%saveas(fig,['dataset' mat2str(length(dir('*.png'))+1) '.png'])
