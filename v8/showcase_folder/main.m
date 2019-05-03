% main script to run function
% showcase version, just for show!

clear; close all; clc; tic;

%% Set parameters 

% System params %
Parameters.tmax     = 3e4;     % maximum time of simulation
Parameters.dt       = 0.2;     % time step

% Neurons params %

Parameters.iu1      = -0.0009; % mean I parameter for first population
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
Parameters.istate   = 1;
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

Options.doAplot       = 1;
Options.doplot1       = 0; Options.doraster = 1; %(0-none,1-yes,2-colors)
Options.doplot2       = 1;
Options.dogifplot     = 0;

%% Call simulation
simulate_v1(Parameters, Options);

toc
