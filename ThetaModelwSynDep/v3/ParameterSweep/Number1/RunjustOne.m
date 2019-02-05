% Run just one simulation with plots

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
isig1 = 0.0;

tmax = 1000;

istate = 2; % state to start in (1 = low, 2 = high)
bumpit = 1; % should we perturb system during simulation? (give it a kick)
doplots = 1;% want to plot anything?

D = [2.71794871794872];
%D = 3;


synctheta_v3PS(D,ydrop,tauy,iu1,isig1,istate,bumpit,doplots,tmax);





























