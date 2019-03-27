% Main script to run other scripts for theta model v5

clear; close all; clc

% Params %

params.n1 = 135;   % number of neurons in the first population
params.n2 = 15;     % number of neurons in the second population
params.dt = 0.05;   % time step
params.tmax = 3e4;     % maximum time of simulation
params.iu1 = -0.15;  % mean I parameter for first population
params.isig1 = 0.000;  % std of I parameter for first population
params.iu2 = 0.01;  % mean I parameter for second population#
params.isig2 = 0.000;    % std of I parameter for second population
params.prob = 0.5; % E-R graph, prob is prob of connection.
params.tauavg=1e2;   % Relaxation of network excitement

params.ydrop = .01; % How much of an affect firing has on synaptic depression
% (should be between 0 and 1)!!!
params.tauy  =  5e3; % Char time for return to ss for y (synap depress)

params.silence=0;
params.bumpit=0;
params.istate=3;

%% First param sweep

if 0
    ParamsweepMain
else
    load PSresults7b
end

%% Get synaptic depression nulcline

if 0
    params.tauy  =  5e3;
    params.ydrop = .01;
    params.MaxD = 3;
    syndepresnulc_main
else
    load SDn_results7b.mat
end

%% Simulate full system (w/ syn dep back in)

if 1
    
    params.tmax = 3e4;

h = figure('Position',[1000 400 1300 900]); subplot(5,3,[1 6]);
plot(Darray/params.MaxD,Lowsps); hold on;
plot(Darray/params.MaxD,Highsps); 
plot(yssarray,spsarray,'bx')
legend('Low init state','High init state + mid simulation bump','y_{ss}');
title('Strength of synapse determines possible fates of network')
xlabel('network-wide mean y (Average network synaptic depression)')
ylabel('Spikes/Second (of network*)')

subplot(5,3,[13 15]); hold on;
str = ['isig1 = ' mat2str(params.isig1) '. iu1 = ' mat2str(params.iu1) ...
    '. isig2 = ' mat2str(params.isig2) '. iu2 = ' mat2str(params.iu2) ...
    '. n1 = ' mat2str(params.n1) '. n2 = ' mat2str(params.n2) '. prob = .5' ...
    '. tauy = ' mat2str(params.tauy) '. ydrop = ' mat2str(params.ydrop)];
annotation('textbox',[.2 .08 .1 .1],'String',str,'FitBoxToText','on');

str = ['tmax = ' mat2str(params.tmax) '. DMax = ' mat2str(params.MaxD)];
annotation('textbox',[.2 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

synctheta_v5(params,h)

end