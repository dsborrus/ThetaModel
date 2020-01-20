% main script for creating and/or calling Figure Object
% this directory is for testing network topoligies and comparing them to
% prebotc


% for bioloigical equivelent experiments see
% Hayes et al 2012 in PNAS
% Wang et al 2014 in eLife
clear; close all; clc; t3t = tic;
% = 5/30/19

disp('Running main script for testing network topoligies.');

%% Params (Note, these params only change the initial networks. After they 
% are built, you need to change params inside the methods function of FigureArms

% conmat params
Parameters.tmax     = 100000;     % maximum time of simulation (in ms)
Options.conmat      = 0;
    % case 1 - ER
Parameters.er_prob  = 0.2;     % prob of connection
    % case 2 - small world
Parameters.sw_M     = 7;       % small world - number of Ns on each side
Parameters.sw_p     = 0.05;     % small world - probability of "short cut" 
    % case 3 - scale-free
Parameters.sf_mo    = 30;      % scale free - size of seed
Parameters.sf_m     = 30;      % sacle free - average degree (use mo=m or m<mo)
Parameters.sf_d     = 1;     % d = chance that a wiring will be reciprical
Options.startstruct.s= 0;      % starting structure (0 - clique)(1 - ER)
Options.startstruct.p=.6;     % p of ER network, if making ER 
    % case 4 - premade scale free
%Parameters.A2load = 'saved600neuronscalefree.mat';
    % case 5 - Klemm and Eguilez ( 2002b) scalefree, small world
Parameters.ke_mo = 25;          % scale free - size of seed
Parameters.ke_mu = .5;         % The chance that a wiring will wire to a non active node
Parameters.ke_d = 0.0;         % d = chance that a wiring will NOT be reciprical
    

% ablation params

% ablation params %

Options.Doablate    = 1;
Parameters.t2ablat  = 250;     % When to ablate neurons (every XXX seconds)
Parameters.N2k_w1p  = 5;       % number of neurons to kill with one pulse

% Other parameters 

% System params %
Parameters.dt       = 0.25;     % time step

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

Options.doAplot       = 0;
Options.doplot1       = 0; Options.doraster = 1; %(0-none,1-yes,2-colors)
Options.doplot2       = 0;
Options.dogifplot     = 0;
Options.trackstatevariablesmeans = 1;

%% Body

% str1 = 'ER_p20';
% str2 = 'sw_M7_p10';
% str3 = 'sw_M7_p40';
% str4 = 'sf_mmo30_d100_ssclique';
% str5 = 'sf_m40_mo20_d030_ssERp050';
% str6 = 'ke_mo25_mu020_d000_ssclique';
% str7 = 'ke_mo50_mu040_d080_ssERp050';

counter = 0;    % keeps track of number of networks

%% Make or Load Networks

disp('making or loading networks')

% ER w/ prob 014
if 0
    str1 = 'ER_p14';
    Options.conmat      = 1;
    Parameters.er_prob  = 0.14;
    F1 = FigureArms(Parameters,Options,str1);
    F1.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F1.mat','F1')
else
    load('BackUpNetworks/F1.mat')
end

% small world with M 7 and p 0.1
if 0
    str2 = 'sw_M7_p10';
    Options.conmat      = 2;
    Parameters.sw_M     = 7;
    Parameters.sw_p     = 0.10;
    F2 = FigureArms(Parameters,Options,str2); 
    F2.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F2.mat','F2')
else
    load('BackUpNetworks/F2.mat','F2')
end

% small world with M 7 and p 0.4
if 0
    str3 = 'sw_M7_p40';
    Options.conmat      = 2;
    Parameters.sw_M     = 7;
    Parameters.sw_p     = 0.40;
    F3 = FigureArms(Parameters,Options,str3);
    F3.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F3.mat','F3')
else
    load('BackUpNetworks/F3.mat','F3')
end

% scale free with M & Mo 30, d=1, and start is clique
if 0
    str4 = 'sf_mmo7_d100_ssclique';
    Options.conmat      = 3;
    Parameters.sf_mo    = 7; 
    Parameters.sf_m     = 7;    
    Parameters.sf_d     = 1.00;
    Options.startstruct.s= 0;
    Options.startstruct.p=.6;
    F4 = FigureArms(Parameters,Options,str4);
    F4.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F4.mat','F4')
else
    load('BackUpNetworks/F4.mat','F4')
end

% scale free with M 7 & Mo 14, d=0.0, and start is ER with p.5
if 0
    str5 = 'sf_m7_mo14_d000_ssERp050';
    Options.conmat      = 3;
    Parameters.sf_mo    = 14; 
    Parameters.sf_m     = 7;    
    Parameters.sf_d     = 0.0;
    Options.startstruct.s= 1;
    Options.startstruct.p=.50;
    F5 = FigureArms(Parameters,Options,str5);
    F5.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F5.mat','F5')
else
    load('BackUpNetworks/F5.mat','F5')
end

% scale-free, small-world with Mo7, mu=.50, d=0.00, and start is ER prob
% 0.14
if 0
    str6 = 'ke_mo7_mu050_d000_ssERp014';
    Options.conmat   = 5;
    Parameters.ke_mo = 7;          % scale free - size of seed
    Parameters.ke_mu = .50;         % The chance that a wiring will wire to a non active node
    Parameters.ke_d = 0.0; 
    Options.startstruct.s= 1;
    Options.startstruct.p=0.14;
    F6 = FigureArms(Parameters,Options,str6);
    F6.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F6.mat','F6')
else
    load('BackUpNetworks/F6.mat','F6')
end

% scale-free small-world with Mo50, mu=.40, d=0.80, and start is ER p.5
if 0
    str7 = 'ke_mo50_mu040_d080_ssERp050';
    Options.conmat   = 5;
    Parameters.ke_mo = 50;          % scale free - size of seed
    Parameters.ke_mu = .40;         % The chance that a wiring will wire to a non active node
    Parameters.ke_d = 0.80; 
    Options.startstruct.s= 1;
    Options.startstruct.p=.5;
    F7 = FigureArms(Parameters,Options,str7);
    F7.MakeGraph;
    counter = counter+1;
    save('BackUpNetworks/F7.mat','F7')
else
    %load('BackUpNetworks/F7.mat','F7')
end

disp('Done making or loading networks')

%% Analyse new or loaded networks

if 0
    disp('Analysing networks')
    F1.AnalyseGraph(1,1);
    F2.AnalyseGraph(0,0);
    F3.AnalyseGraph(0,0);
    F4.AnalyseGraph(2,2);
    F5.AnalyseGraph(2,2);
    F6.AnalyseGraph(2,2);
    %F7.AnalyseGraph(2,2);
    disp('Done analysing networks')
else
    disp('Skipping network connectivity analysis')
end

%% Simulate networks

if 0

    disp('Simulating networks, without ablation')

    if 1
        loadlastsim = 1;
        F1.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F1.mat','F1')
        disp('Simmed network 1')
    end

    if 1
        loadlastsim = 1;
        F2.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F2.mat','F2')
        disp('Simmed network 2')
    end

    if 1
        loadlastsim = 1;
        F3.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F3.mat','F3')
        disp('Simmed network 3')
    end

    if 1
        loadlastsim = 1;
        F4.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F4.mat','F4')
        disp('Simmed network 4')
    end

    if 1
        loadlastsim = 1;
        F5.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F5.mat','F5')
        disp('Simmed network 5')
    end

    if 1
        loadlastsim = 1;
        F6.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F6.mat','F6')
        disp('Simmed network 6')
    end

    if 0
        loadlastsim = 0;
        F7.simplsimNetwork(loadlastsim)
        save('BackUpNetworks/F7.mat','F7')
    end

    disp('Networks have been simulated simply (or have been loaded from saves) and plotted')
else
    disp('Didnt run simple,intact,unablated sims. They must already be loaded.')
end

%% Simulate and ablate network
if 0

    disp('Simulating networks, WITH ablation')

    % loadlaststim allows us to skip the simulation part if we only care
    % about the graphic output. But if we changed settings/parameters, we
    % should rerun the sim!!
    if 1
        loadlastsim = 1;
        F1.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F1.mat','F1')
        disp('Ablated network 1')
    end

    if 1
        loadlastsim = 1;
        F2.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F2.mat','F2')
        disp('Ablated network 2')
    end

    if 1
        loadlastsim = 1;
        F3.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F3.mat','F3')
        disp('Ablated network 3')
    end

    if 1
        loadlastsim = 1;
        F4.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F4.mat','F4')
        disp('Ablated network 4')
    end

    if 1
        loadlastsim = 1;
        F5.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F5.mat','F5')
        disp('Ablated network 5')
    end

    if 1
        loadlastsim = 1;
        F6.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F6.mat','F6')
        disp('Ablated network 6')
    end

    if 0
        loadlastsim = 0;
        F7.simublateNetwork(loadlastsim)
        save('BackUpNetworks/F7.mat','F7')
        disp('Ablated network 7')
    end

    disp('Networks have been ablated (or have been loaded from saves) and plotted')
else
    disp('Didnt run ablation sims. They must already be loaded.')
end

%% Zoom in and plot some ablation traces

if 0
    disp('Zooming in and plotting ablation traces')
    window = 60;
    windowfrontbuffer = 200;
    
    for chunknumber = [1 2 3 4]
    F1.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    F2.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    F3.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    F4.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    F5.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    F6.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    %F7.plotchunksofablation(chunknumber,window,windowfrontbuffer)
    end
    
end

%% Rhythm summary statistics

if 0
   
    disp('Comapring rhythm breakdown')
    
    Fs = [F1;F2;F3;F4;F5;F6];
    
    F1.AblationRhythmSummary(Fs);
    
end

%% latex stuff

fileID = fopen('tex_v2/defs.tex','w');
fprintf (fileID,['\\def \\nNetworks{' num2str(counter) '}\n']);
fprintf (fileID,['\\def \\strone{' F1.str '}\n']);
fprintf (fileID,['\\def \\strtwo{' F2.str '}\n']);
fprintf (fileID,['\\def \\strthree{' F3.str '}\n']);
fprintf (fileID,['\\def \\strfour{' F4.str '}\n']);
fprintf (fileID,['\\def \\strfive{' F5.str '}\n']);
fprintf (fileID,['\\def \\strsix{' F6.str '}\n']);
%fprintf (fileID,['\\def \\strseven{' F7.str '}\n']);
fprintf (fileID,['\\def \\neuronkill{' mat2str(F1.Parameters.N2k_w1p) '}\n']);

% system('pdflatex tex_v2/MainTex.tex')

%% State variable dynamics

if 0
    disp('State variable dynamics')
    
    if 1
    % make sure chunks is a row vector!!
    chunks = [0 3];
    vars2comp = {'f','s'};
    figindex = 1;
    
    F1.DrawStateTraces(chunks,vars2comp,figindex)
    F2.DrawStateTraces(chunks,vars2comp,figindex)
    F3.DrawStateTraces(chunks,vars2comp,figindex)
    F4.DrawStateTraces(chunks,vars2comp,figindex)
    F5.DrawStateTraces(chunks,vars2comp,figindex)
    F6.DrawStateTraces(chunks,vars2comp,figindex)
    end
    
    if 1
    chunks = [0 3];
    vars2comp = {'y','s'};
    figindex = 2;
    
    F1.DrawStateTraces(chunks,vars2comp,figindex)
    F2.DrawStateTraces(chunks,vars2comp,figindex)
    F3.DrawStateTraces(chunks,vars2comp,figindex)
    F4.DrawStateTraces(chunks,vars2comp,figindex)
    F5.DrawStateTraces(chunks,vars2comp,figindex)
    F6.DrawStateTraces(chunks,vars2comp,figindex)
    end
    
    if 1
    chunks = [0 3];
    vars2comp = {'y','f'};
    figindex = 3;
    
    F1.DrawStateTraces(chunks,vars2comp,figindex)
    F2.DrawStateTraces(chunks,vars2comp,figindex)
    F3.DrawStateTraces(chunks,vars2comp,figindex)
    F4.DrawStateTraces(chunks,vars2comp,figindex)
    F5.DrawStateTraces(chunks,vars2comp,figindex)
    F6.DrawStateTraces(chunks,vars2comp,figindex)
    end 
end

%% Network Statistics during ablation

if 1

    disp('Calculating network statistics')
    % I believe we want a row vector here
    Fs = [F1,F2,F3,F4,F5,F6];
    F1.NetworkStatisticsDuringAblation(Fs);
    
    
    
end