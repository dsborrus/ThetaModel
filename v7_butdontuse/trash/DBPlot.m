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


ax2 = subplot(4,1,2); hold on;
for i = 1:n
    plot(t(tw),raster((tw),i)*i,'.k'); hold on;
end


%mr = max(max(raster));
%map = [linspace(.5,1,mr)'  linspace(.5,1,mr)' linspace(.5,1,mr)'];
% workinraster = raster;   
% 
% temp = 0;
% while any(any(workinraster)) ~= 0   
%     temp = temp+1;
%     for i = 1:n
%        plot(t(tw),workinraster((tw),i)*i,'.','color',map(temp,:)); hold on;
%     end
%     workinraster = workinraster - 1;
%     workinraster(workinraster<0) = 0;
% end


title('Raster Plot')

ax3 = subplot(4,1,3);
rr = randi(n1);
plot(t(tw),y(tw,rr))
title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3],'x')

subplot(4,1,4)
hold on
str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off
