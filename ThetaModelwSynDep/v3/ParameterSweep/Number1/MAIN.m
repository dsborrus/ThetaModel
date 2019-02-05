% Script to run parameter sweep

% Version 2 - removing log stuff, going linear


% Version 1 - log plotting

% First Parameters sweep

% Goal is to find the frequency of the network as a function of D (strength
% of pulse), then generate a nice fig from that

% Expecting that at low D, network does not sustain activity
% at high D, network is in completely active state
% at medium D, there is bistability between down and up state

close all;
clear;
clc
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);
tic

ydrop = 0;
tauy = 8;
iu1 = -0.01;
isig1 = 0.007;

tmax = 1000;

%istate = 2; % state to start in (1 = low, 2 = high)
%bumpit = 1; % should we perturb system during simulation? (give it a kick)
doplots = 0;% want to plot anything?

% Ok here's the plan. Going to try a range of Ds. with tmax = 1e3? Will run
% each D value twice. Once with low initial state and no bump. and again
% with high initial state and bump. See if we can find a bistable range of
% values.
Dstart = 1;
Dend = 5;
nD = 200;
D = linspace(Dstart,Dend,nD);

flowout = zeros(nD,1);
fhighout = zeros(nD,1);

for d = 1:nD

    istate = 1;
    bumpit = 0;
    flowout(d) = synctheta_v3PS(D(d),ydrop,tauy,iu1,isig1,istate,bumpit,doplots,tmax);

    istate = 2;
    bumpit = 1;
    fhighout(d) = synctheta_v3PS(D(d),ydrop,tauy,iu1,isig1,istate,bumpit,doplots,tmax);
    
    disp(['Done ' mat2str(d) ' out of ' mat2str(nD)]);

end

%% figure stuff

figure('Position',[1200 500 900 700]); subplot(3,3,[1 6]);
plot(D,flowout); hold on;
plot(D,fhighout); legend('Low init state','High init state + mid simulation bump');
title('Strength of synapse determines possible fates of network')
xlabel('D')
ylabel('Spikes/Second (of network*)')

subplot(3,3,[7 9]); hold on;
str = ['isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. n = 200. prob = .5' ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.2 .2 .1 .1],'String',str,'FitBoxToText','on');

axis off




























