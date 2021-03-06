% Main script to run other simulation - ablation studies

% main script to run function
% showcase version, just for show!

clear; close all; clc; tic; rng(1000);

addpath('..')

%% conmat params
Parameters.tmax     = 3e4;     % maximum time of simulation
Options.conmat      = 1;
    % case 1 - ER
Parameters.er_prob  = 0.2;     % prob of connection
    % case 2 - small world
Parameters.sw_M     = 7;       % small world - number of Ns on each side
Parameters.sw_p     = 0.05;     % small world - probability of "short cut" 
    % case 3 - scale-free
Parameters.sf_mo    = 30;      % scale free - size of seed
Parameters.sf_m     = 30;      % sacle free - average degree (use mo=m or m<mo)
Parameters.sf_d     = 1.30;     % d = chance that a wiring will be reciprical
Options.startstruct.s= 0;      % starting structure (0 - clique)(1 - ER)
Options.startstruct.p=.6;     % p of ER network, if making ER 
    % case 4 - premade scale free
Parameters.A2load = 'saved600neuronscalefree.mat';
    % case 5 - Klemm and Eguilez ( 2002b) scalefree, small world
Parameters.ke_mo = 25;          % scale free - size of seed
Parameters.ke_mu = .5;         % The chance that a wiring will wire to a non active node
Parameters.ke_d = 0.0;         % d = chance that a wiring will NOT be reciprical
    

%% ablation params

% ablation params %

Options.Doablate    = 0;
Parameters.t2ablat  = 10;      % When to ablate neurons (every XXX seconds)
Parameters.N2k_w1p  = 5;       % number of neurons to kill with one pulse

%% Other parameters 

% System params %
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

% Script params %

Options.doAplot       = 1;
Options.doplot1       = 0; Options.doraster = 1; %(0-none,1-yes,2-colors)
Options.doplot2       = 1;
Options.dogifplot     = 0;

%% Call simulation
Output = simulate_v1(Parameters, Options);

toc
