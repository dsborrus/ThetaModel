function DBPlot_v2(dt,tmax,t,y,m,n,sij,spikes,raster,rr,Ihistory,sijhistory,vin,istate)

D = vin(1);
isig1 = vin(2);
iu1 = vin(3);
isig2 = vin(4);
iu2 = vin(5);
N = vin(6);
prob = vin(7);
tauavg = vin(8);
mgain = vin(9);
taum = vin(10);
nrise = vin(11);
taun = vin(12);
sigain = vin(13);
tausi = vin(14);
tautheta=vin(15);
noisesigma=vin(16);
conmat = vin(17);
doablate = vin(18);
t2ablat = vin(19);
numN2kill_wonepulse = vin(20);

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

nplots = 9;

tstart = dt;
tend = tmax;

tw = tstart/dt:tend/dt;

fig = figure('Position',[800 500 1500 1200]);

ax1 = subplot(nplots,1,1);
hold on
plot(t(tw)/1000,(spikes(tw)./tauavg)*1000)
plot(t2ablat:t2ablat:tmax/1000,zeros(length(t2ablat:t2ablat:tmax/1000),1)+10,'ro')
title('Network Activity')


ax2 = subplot(nplots,1,2); hold on;
for i = 1:N
    plot(t(tw)/1000,raster((tw),i)*i,'.k'); hold on;
end

title('Raster Plot')

if 1
ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,mean(y(tw,:),2))
title(['Mean synaptic depression (y) (n=' mat2str(rr) ')'])
else
ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,y(tw,rr))
title(['One random neuron`s synaptic depression (y) (n=' mat2str(rr) ')'])
end

ax4 = subplot(nplots,1,4);
plot(t(tw)/1000,sij(tw,rr))
title(['One random neuron`s synaptic current output (s) (n=' mat2str(rr) ')'])

ax5 = subplot(nplots,1,5);
plot(t(tw)/1000,m(tw,rr))
title(['One random neuron`s (m) (n=' mat2str(rr) ')'])

ax6 = subplot(nplots,1,6);
plot(t(tw)/1000,n(tw,rr))
title(['One random neuron`s (n) (n=' mat2str(rr) ')'])

ax7 = subplot(nplots,1,7);
plot(t(tw)/1000,sijhistory(tw))
title(['One random neuron`s synaptic current input (n=' mat2str(rr) ')'])

ax8 = subplot(nplots,1,8);
plot(t(tw)/1000,Ihistory(tw))
title(['One random neuron`s total applied current (I+delta*syncurrent*syndepression) (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'x')

subplot(nplots,1,nplots);
hold on

str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. N = ' mat2str(N) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. mgain = ' mat2str(mgain) '. taum = ' mat2str(taum) ...
    '. sigain = ' mat2str(sigain) '. tausi = ' mat2str(tausi)];
annotation('textbox',[.1 .08 .1 .1],'String',str,'FitBoxToText','on');

str = ['nrise = ' mat2str(nrise) '. taun = ' mat2str(taun) ...
       '. tautheta = ' mat2str(tautheta) '. noisesigma = ' mat2str(noisesigma)...
       '. istate = ' mat2str(istate) '. conmat = ' mat2str(conmat) '. doablate = ' mat2str(doablate)...
       '. t2ablat = ' mat2str(t2ablat) '. numN2kill_wonepulse = ' mat2str(numN2kill_wonepulse)];
annotation('textbox',[.1 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

pngindx = length(dir('*.png'))+1;

saveas(fig,['dataset' mat2str(pngindx) '.png'])

save(['dataset' mat2str(pngindx) '.mat'],'spikes')
