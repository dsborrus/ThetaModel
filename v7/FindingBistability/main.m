% Main script to run other scripts for theta model v5

clear; close all; clc

% Params %

params.n1 = 100;   % number of neurons in the first population
params.n2 = 00;     % number of neurons in the second population
params.dt = 0.5;   % time step
params.tmax = 1e4;     % maximum time of simulation
params.iu1 = 0;  % mean I parameter for first population
params.isig1 = 0.0;  % std of I parameter for first population
params.iu2 = 1e-6;  % mean I parameter for second population#
params.isig2 = 0.000;    % std of I parameter for second population
params.prob = 1.0; % E-R graph, prob is prob of connection.
params.tauavg = 1e2;   % Relaxation of network excitement

params.sigain = .5;
params.tausi = 5;

%params.ydrop = .01; % How much of an affect firing has on synaptic depression
% (should be between 0 and 1)!!!
%params.tauy  =  5e3; % Char time for return to ss for y (synap depress)

%params.silence=0;
%params.bumpit=0;
%params.istate=3;

% %% Start cluster
% poolid = gcp('nocreate');
% if isempty(poolid)
% c = parcluster;
% c.NumWorkers=8;
% parpool(c);
% end

%% First param sweep

if 1
    Main_Paramsweep
else
    load PSresults5
end

%% Get synaptic depression nulcline

if 1
    params.tauy  =  5e3;
    params.ydrop = .1;
    params.MaxD = .061;
    Main_syndepresnulc
else
    load SDn_results5
end

%% Simulate full system (w/ syn dep back in)

if 1
    params.silence=0;
    params.bumpit=0;
    params.istate=3;
    params.tmax = 2e4;

h = figure('Position',[1000 400 1300 900]); subplot(5,3,[1 6]);
plot(Darray/params.MaxD,Lowsps,'linewidth',2); hold on;
plot(Darray/params.MaxD,Highsps,'linewidth',2); 
plot(yssarray,spsarray,'--','color',[0,0,1]);
plot(yssn2_array,spsarray,'--','color',[0,0.5,0]);
legend('Low init state','High init state + mid simulation bump','y_{ss}');
title('Strength of synapse determines possible fates of network')
xlabel('network-wide mean y (Average network synaptic depression)')
ylabel('Spikes/Second (of network*)')

subplot(5,3,[13 15]); hold on;
str = ['isig1 = ' mat2str(params.isig1) '. iu1 = ' mat2str(params.iu1) ...
    '. isig2 = ' mat2str(params.isig2) '. iu2 = ' mat2str(params.iu2) ...
    '. n1 = ' mat2str(params.n1) '. n2 = ' mat2str(params.n2) '. prob = ' mat2str(params.prob) ...
    '. tauy = ' mat2str(params.tauy) '. ydrop = ' mat2str(params.ydrop)];
annotation('textbox',[.2 .08 .1 .1],'String',str,'FitBoxToText','on');

str = ['tmax = ' mat2str(params.tmax) '. DMax = ' mat2str(params.MaxD)];
annotation('textbox',[.2 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

synctheta_v6(params,h)

end