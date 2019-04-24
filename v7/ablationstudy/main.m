% Main script to run other scripts for theta model v5

clear; close all; clc

tic

% A is connectivity matrix (1-ER, 2-small world, 3-scale free, 4-hierarchy)
A = 1;
ablate = 0;

synctheta_v7_edit(2e4,1,1,0,A,ablate);

toc
