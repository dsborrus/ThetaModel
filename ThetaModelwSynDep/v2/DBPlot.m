% DB plot

figure

ax1=subplot(4,1,1);
hold on
plot(t,spikes)

ax2=subplot(4,1,2);
hold on
plot(t,thetaIMP(:,1))

ax3=subplot(4,1,3);
hold on
plot(t,yIMP(:,1))

linkaxes([ax1, ax2, ax3],'x')

subplot(4,1,4)
hold on
str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off


% Ron wuz here