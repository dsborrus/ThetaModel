% Run just one simulation with plots

% good if you want to see the actual output of a simulation

close all;
clear;
clc
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);
tic

% % % % % % % %  Params % % % % % % % % 
ydrop = 0;
tauy = 8;
n1 = 540;
n2 = 60;
iu1 = -0.01;
isig1 = 0.000;
iu2 = 0.01;
isig2 = 0.000;

tmax = 2000;

% WHICH STATE TO START IN?
istate = 2; % state to start in (1 = low, 2 = high)
% SHOULD WE PERTURB SYSTEM?
bumpit = 1; % should we perturb system during simulation? (give it a kick)

% plotting? Probably.
doplots = 1;% want to plot anything?

% % % %  important param % % % % 
D = 3;
%D = 3;

% % % %  RUN % % % % % % % % 

synctheta_v4PS(D,ydrop,tauy,n1,n2,iu1,isig1,iu2,isig2,istate,bumpit,doplots,tmax);

toc
























