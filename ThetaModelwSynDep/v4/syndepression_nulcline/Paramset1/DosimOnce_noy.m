% script to sim
% network activity, without synaptic depression

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

ydrop = .1;
tauy = 15;

tmax = 2e3; % maximum time for simulation. Below 1000 you might start 
            % analysing too early (we analyse last half of the simulation)
            
plotit = 1;
            
D = 3;


istate = 1;
bumpit = 0;
synctheta_function_noy(D,ydrop,tmax,tauy,istate,bumpit,n1,n2,iu1,iu2,isig1,isig2,plotit);
    
