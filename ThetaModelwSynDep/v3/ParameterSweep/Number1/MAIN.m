% Script to Run: parameter sweep

% % % % Version log % % % % 

% Version 2 - removing log stuff, going linear

% Version 1 - log plotting

% % % % %  Description % % % %  % % % % 
% First Parameters sweep

% Goal is to find the frequency of the network as a function of D (strength
% of pulse), then generate a nice fig from that

% Expecting that at low D, network does not sustain activity
% at high D, network is in completely active state
% at medium D, there is bistability between down and up state

% note for user. If you want to see a single simulation. Run "RunjustOne.m"
% or go back to /v3 and run synchtheta_v3C

% % % % Random pre-stuff % % % % % % % % 

close all;
clear;
clc
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);
tic

% % % % Params % % % % % % % % 

ydrop = 0;  % set synaptic depression to 0. Dont wan't D changing "slow"ly during simulation. Or do you??
tauy = 8;
n1 = 180;
n2 = 20;    % 10% of population has intrinsic rhythm.
iu1 = -0.01;
isig1 = 0.000;
iu2 = 0.01; % See v3/I2freq for mapping between I and freq
isig2 = 0;

tmax = 1000; % maximum time for simulation. Below 1000 you might start 
             % analysing to early (we analyse last half of the simulation)

doplots = 0; % want to plot simulations? Probably not for a param sweep

% % % %  Ok here's the good stuff.  % % % %  % % % % 

% Going to try a range of Ds. with tmax = 1e3. Will run
% each D value twice. Once with low initial state and no bump. and again
% with high initial state and bumps. See if we can find a bistable range of
% values.
Dstart = 1;
Dend = 4;
nD = 200;
D = linspace(Dstart,Dend,nD);

% % % % Script really starts here % % % % 

% initialize output
flowout = zeros(nD,1);
fhighout = zeros(nD,1);

for d = 1:nD

    % sim down state
    istate = 1;
    bumpit = 0;
    flowout(d) = synctheta_v4PS(D(d),ydrop,tauy,n1,n2,iu1,isig1,iu2,isig2,istate,bumpit,doplots,tmax);

    % sim up state
    istate = 2;
    bumpit = 1;
    fhighout(d) = synctheta_v4PS(D(d),ydrop,tauy,n1,n2,iu1,isig1,iu2,isig2,istate,bumpit,doplots,tmax);
    
    disp(['Done ' mat2str(d) ' out of ' mat2str(nD)]);

end

%% figure stuff

figure('Position',[1200 400 1100 700]); subplot(3,3,[1 6]);
plot(D,flowout); hold on;
plot(D,fhighout); legend('Low init state','High init state + mid simulation bump');
title('Strength of synapse determines possible fates of network')
xlabel('D')
ylabel('Spikes/Second (of network*)')

subplot(3,3,[7 9]); hold on;
str = ['isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. n1 = ' mat2str(n1) '. n2 = ' mat2str(n2) '. prob = .5' ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.2 .2 .1 .1],'String',str,'FitBoxToText','on');

axis off




























