% DB plot %
close all

tstart = 5000;

tend = 5500;

if tstart>=tmax
    tstart = dt;
end

if tend>tmax
    tend = tmax;
end

tw = tstart/dt:tend/dt;

figure('Position',[800 500 1300 800])

ax1 = subplot(4,1,1);
hold on
plot(t(tw),spikes(tw))
title('Network Activity')


ax2 = subplot(4,1,2);
temp = dt:dt:tmax;
for i = 1:n
    plot(temp(tw),raster((tw),i)*i,'k.'); hold on;
end
title('Raster Plot')

ax3 = subplot(4,1,3);
rr = randi(n);
plot(t(tw),y(tw,rr))
title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3],'x')

subplot(4,1,4)
hold on
str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off
