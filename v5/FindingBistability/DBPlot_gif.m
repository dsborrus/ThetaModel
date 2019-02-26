function DBPlot_v2(dt,tmax,t,n,y,sij,spikes,raster,rr,Ihistory,sijhistory,vin)

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

figure('position',[500 500 1000 600])
subplot(6,3,[1 6])

%% figs %%

temp = length(dir('*.gif'));
filename = (['AnimatedGIF' mat2str(temp) '.gif']);

h = figure('Position',[1000 400 1300 900]); subplot(5,3,[1 6]);
% plot(D,flowout); hold on;
% plot(D,fhighout); 
%plot(yss_array*DMax,sps_array,'bx')
plot(D/DMax,flowout); hold on;
plot(D/DMax,fhighout); 
plot(yss_array,sps_array,'bx')
legend('Low init state','High init state + mid simulation bump','y_{ss}');
title('Strength of synapse determines possible fates of network')
xlabel('network-wide mean y (Average network synaptic depression)')
ylabel('Spikes/Second (of network*)')

subplot(5,3,[13 15]); hold on;
str = ['isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. n1 = ' mat2str(n1) '. n2 = ' mat2str(n2) '. prob = .5' ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.2 .08 .1 .1],'String',str,'FitBoxToText','on');

str = ['tmax = ' mat2str(tmax) '. DMax = ' mat2str(DMax)];
annotation('textbox',[.2 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

subplot(5,3,[7 8]);
pl1 = plot(0,0);
pl2 = plot(0,0);
pl3 = plot(0,0);
ylim([0 6e-3]);
axis off
subplot(5,3,[10 11]);

wind = 30; %seconds
pause(1)
start = tmax-105;
for i = start:tmax-5
    subplot(5,3,[1 6]); hold on;
    delete(pl1);
    delete(pl3)
%     pl1 = plot(yss_tarray(i)*DMax,sps_tarray(i),'ro');
%     pl3 = plot(yss_tarray(i-2:i)*DMax,sps_tarray(i-2:i),'r');
     pl1 = plot(yss_tarray(i),sps_tarray(i),'ro');
     pl3 = plot(yss_tarray(i-2:i),sps_tarray(i-2:i),'r');

        
    subplot(5,3,[7 8]); hold on;
    delete(pl2)
    pl2 = plot(t((i-wind)/dt:(i/dt)),spikes((i-wind)/dt:(i/dt)),'b');
    
    subplot(5,3,[10 11]); hold off;
    for k = 1:(n1+n2)
        plot(t((i-wind)/dt:(i/dt)),raster(((i-wind)/dt:(i/dt)),k)*k,'k.'); hold on;
    end
    
    drawnow
    %pause(1)
    
        if dogif
        % capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write the GIF File
        if i == start
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append')
        end
        end
    
end

% nplots = 7;
% 
% tstart = dt;
% tend = tmax;
% 
% tw = tstart/dt:tend/dt;
% 
% figure('Position',[800 500 1300 800])
% 
% ax1 = subplot(nplots,1,1);
% hold on
% plot(t(tw)/1000,spikes(tw))
% title('Network Activity')
% 
% 
% ax2 = subplot(nplots,1,2); hold on;
% for i = 1:n
%     plot(t(tw)/1000,raster((tw),i)*i,'.k'); hold on;
% end
% 
% title('Raster Plot')
% 
% ax3 = subplot(nplots,1,3);
% plot(t(tw)/1000,y(tw,rr))
% title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])
% 
% ax4 = subplot(nplots,1,4);
% plot(t(tw)/1000,sij(tw,rr))
% title(['One random neuron`s synaptic current output (n=' mat2str(rr) ')'])
% 
% ax5 = subplot(nplots,1,5);
% plot(t(tw)/1000,sijhistory(tw))
% title(['One random neuron`s synaptic current input (n=' mat2str(rr) ')'])
% 
% ax6 = subplot(nplots,1,6);
% plot(t(tw)/1000,Ihistory(tw))
% title(['One random neuron`s total applied current (I+delta*syncurrent*syndepression) (n=' mat2str(rr) ')'])
% 
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
% 
% subplot(nplots,1,nplots)
% hold on
% 
% D = vin(1);
% isig1 = vin(2);
% iu1 = vin(3);
% isig2 = vin(4);
% iu2 = vin(5);
% n = vin(6);
% prob = vin(7);
% tauavg = vin(8);
% tauy = vin(9);
% ydrop = vin(10);
% 
% str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
%     '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
%     '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
%     '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
% annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');
% 
% axis off
end