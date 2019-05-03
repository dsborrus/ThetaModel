% Main script to run other scripts for theta model v5

clear; close all; clc

tic

[spikes,dt] = synctheta_v7(1e5,-.00097,4,0,0);

[meanperiod, meanAMP, brstbrstltratio] = psanalysis(spikes,dt,1e5);

toc
