function DBPlot_v2(dt,tmax,t,n,y,sij,spikes,raster,rr,Ihistory,sijhistory,vin)

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

nplots = 7;

tstart = dt;
tend = tmax;

tw = tstart/dt:tend/dt;

figure('Position',[800 500 1300 800])

ax1 = subplot(nplots,1,1);
hold on
plot(t(tw)/1000,spikes(tw))
title('Network Activity')


ax2 = subplot(nplots,1,2); hold on;
for i = 1:n
    plot(t(tw)/1000,raster((tw),i)*i,'.k'); hold on;
end

title('Raster Plot')

ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,y(tw,rr))
title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])

ax4 = subplot(nplots,1,4);
plot(t(tw)/1000,sij(tw,rr))
title(['One random neuron`s synaptic current output (n=' mat2str(rr) ')'])

ax5 = subplot(nplots,1,5);
plot(t(tw)/1000,sijhistory(tw))
title(['One random neuron`s synaptic current input (n=' mat2str(rr) ')'])

ax6 = subplot(nplots,1,6);
plot(t(tw)/1000,Ihistory(tw))
title(['One random neuron`s total applied current (I+delta*syncurrent*syndepression) (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')

subplot(nplots,1,nplots)
hold on

D = vin(1);
isig1 = vin(2);
iu1 = vin(3);
isig2 = vin(4);
iu2 = vin(5);
n = vin(6);
prob = vin(7);
tauavg = vin(8);
tauy = vin(9);
ydrop = vin(10);

str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off
