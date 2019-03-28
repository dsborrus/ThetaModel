% script to run through and create a y_ss line using a shifting D for
% network activity, get bifurcation param w nulcline

close all;
clear;
clc
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);
tic

% % % % Params % % % % % % % % 

load HighResData_Results8.mat

ydrop = .04;
tauy = 28;

tmax = 2e3; % maximum time for simulation. Below 1000 you might start 
            % analysing too early (we analyse last half of the simulation)
            
resolution = 25;            
             
plotit=0;

% DMax would be the parameter to set D, if s = Dmax*y was the pulse.
% However right now it just changes hte position of the y-nulcline in the
% bifurcation figure
DMax = 2.95;
% range of Ds             
% loaded in Dstart and Dend
nD = resolution;
D_new = linspace(Dstart,Dend,nD);

% % % % Script really starts here % % % % 

% initialize output
yss_array = zeros(nD,1);
sps_array = zeros(nD,1);

for d = 1:nD

    % sim up state and record y_ss
    istate = 2;
    bumpit = 1;
    [yss_array(d),sps_array(d)] = synctheta_function_noy(D_new(d),ydrop,tmax,tauy,istate,bumpit,n1,n2,iu1,iu2,isig1,isig2,plotit);
    
    disp(['Done ' mat2str(d) ' out of ' mat2str(nD)]);

end

%% figure stuff

figure
%plot(sps_array,yss_array,'bx')
plot(yss_array,sps_array,'bx')
xlabel('y_ss');ylabel('spikes/second')


figure('Position',[1200 400 1100 700]); subplot(3,3,[1 6]);
plot(D,flowout); hold on;
plot(D,fhighout); 
plot(yss_array*DMax,sps_array,'bx')
legend('Low init state','High init state + mid simulation bump','y_{ss}');
title('Strength of synapse determines possible fates of network')
xlabel('D')
ylabel('Spikes/Second (of network*)')

subplot(3,3,[7 9]); hold on;
str = ['isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. n1 = ' mat2str(n1) '. n2 = ' mat2str(n2) '. prob = .5' ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.2 .2 .1 .1],'String',str,'FitBoxToText','on');

str = ['tmax = ' mat2str(tmax) '. DMax = ' mat2str(DMax)];
annotation('textbox',[.2 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off

save('phaseplaneNulcline.mat')