% main script for changing excitability study
% this main script sets up the param study

clear; close all; clc

addpath('..')

t1=tic;

%% imporatnt params

iu_min = -.00110;
iu_max = -.00080;

resolution = 8;

Parameters.tmax = 3e4;
Parameters.dt = 0.5;

%% standard (not changing) params

% System params %
%Parameters.tmax     = 3e4;     % maximum time of simulation
%Parameters.dt       = 0.2;     % time step

% Neurons params %

%Parameters.iu1      = -0.0009; % mean I parameter for first population
Parameters.isig1    = 0.0;     % std of I parameter for first population
Parameters.iu2      = 0.0;     % mean I parameter for second population#
Parameters.isig2    = 0.0;     % std of I parameter for second population
Parameters.tautheta = 1;       % relaxation of neuron's theta
    % synaptic depression 1
Parameters.mgain    = 0.2;     % How much of a gain firing has on synaptic depression
Parameters.taum     = 300;     % Char time for return to ss for m (synap depress)
    % synaptic depression 2
Parameters.nrise    = 0.011;   % Rate of rise for synaptic depression based on n
Parameters.taun     = 1300;    % Char time for return to ss for n (synap depress)
    % snynaptic conductance
Parameters.sigain   = 1;       % How much of a gain firing has on synaptic conductance
Parameters.tausi    = 15;      % Char time for return to ss for n (synap depress)
    % noise
Parameters.noisesig = 0.009;   % Variance of noise

% Network params %

Parameters.n1       = 100;     % number of neurons in the first population
Parameters.n2       = 0;       % number of neurons in the second population
Parameters.D        = 0.03;    % Strength of networkness
Parameters.tauavg   = 1e2;     % Relaxation of network excitement
Parameters.istate   = 3;
Options.conmat      = 1;
    % case 1 - ER
Parameters.prob     = 0.8;     % prob of connection
    % case 2 - small world
Parameters.sw_M     = 3;       % number of Ns on each side
Parameters.sw_p     = 0.3;     % probability of "short cut" 
    % case 3 - scale-free
Parameters.sf_mo    = 30;      % size of seed
Parameters.sf_m     = 25;      % average degree (use mo=m or m<mo)
    % case 4 - directed

% ablation params %

Options.Doablate    = 0;
Parameters.t2ablat  = 10;      % When to ablate neurons (every XXX seconds)
Parameters.N2k_w1p  = 5;       % number of neurons to kill with one pulse

% Script params %

Options.doAplot       = 0;
Options.doplot1       = 0; Options.doraster = 1; %(0-none,1-yes,2-colors)
Options.doplot2       = 0;
Options.dogifplot     = 0;

%% post-params

iurange = linspace(iu_min,iu_max,resolution);

mp = zeros(resolution,1);
mA = mp;
mbblr = mp;

for k = 1:resolution
    Parameters.iu1 = iurange(k);
    t2 = tic;
    o = simulate_v1(Parameters,Options);
    
    [mp(k), mA(k), mbblr(k)] = psanalysis(o.spikes,Parameters.dt,Parameters.tmax);
    
    disp([mat2str(k) ' out of ' mat2str(resolution) ' done.'])
    toc(t1)
    disp(['Estimated time remaining is ' mat2str(round((toc(t2)/60)*(resolution-k))) ' minutes.'])
    
end

save AAA_data

figure;
subplot(3,1,1);
plot(iurange,mA); title('amplitude')
subplot(3,1,2);
plot(iurange,mp); title('period')
subplot(3,1,3);
plot(iurange,mbblr); title('burst:burstlett ratio')
xlabel('iu values');


toc
