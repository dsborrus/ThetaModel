% Script to run parameter sweep 

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

ydrop = 0.1;
tauy = 8;
iu1 = -0.01;
isig1 = 0.007;

tmax = 1200;

%istate = 2; % state to start in (1 = low, 2 = high)
%bumpit = 1; % should we perturb system during simulation? (give it a kick)
doplots = 0;% want to plot anything?

% Ok here's the plan. Going to try a range of Ds. with tmax = 1e3? Will run
% each D value twice. Once with low initial state and no bump. and again
% with high initial state and bump. See if we can find a bistable range of
% values.
Dstart = 0.5;
Dend = 2;
nD = 12;
D = logspace(log10(Dstart),log10(Dend),nD);

flowout = zeros(nD,1);
fhighout = zeros(nD,1);

for d = 1:nD

    istate = 1;
    bumpit = 0;
    flowout(d) = synctheta_v3PS(D(d),ydrop,tauy,iu1,isig1,istate,bumpit,doplots,tmax);

    istate = 1;
    bumpit = 0;
    fhighout(d) = synctheta_v3PS(D(d),ydrop,tauy,iu1,isig1,istate,bumpit,doplots,tmax);

end

%% figure stuff

figure
semilogx(D,flowout); hold on;
semilogx(D,fhighout); legend('Low state','Highstate');




























