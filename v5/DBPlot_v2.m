function DBPlot_v2(dt,tmax,t,n,n1,y,spikes,raster,vin)

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

tstart = dt;
tend = tmax;

tw = tstart/dt:tend/dt;

figure('Position',[800 500 1300 800])

ax1 = subplot(4,1,1);
hold on
plot(t(tw),spikes(tw))
title('Network Activity')


ax2 = subplot(4,1,2); hold on;
for i = 1:n
    plot(t(tw),raster((tw),i)*i,'.k'); hold on;
end

title('Raster Plot')

ax3 = subplot(4,1,3);
rr = randi(n1);
plot(t(tw),y(tw,rr))
title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3],'x')

subplot(4,1,4)
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
