% script to sim
% network activity, without synaptic depression

close all;
clear;
clc
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);
set(0,'DefaultLegendAutoUpdate','off')
tic

% % % % Params % % % % % % % % 

load phaseplaneNulcline.mat
dt = .01;

%ydrop = .06;
%tauy = 15;

tmax = 2e3; % maximum time for simulation. Below 1000 you might start 
            % analysing too early (we analyse last half of the simulation)
            
plotit = 1;
            
%Dmax = 2.85;


istate = 3; % random starts
bumpit = 0;
[~,~,yss_tarray,sps_tarray,spikes,raster,t] = synctheta_function(DMax,ydrop,tmax,tauy,istate,bumpit,n1,n2,iu1,iu2,isig1,isig2,plotit);

dogif = 0;
    


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