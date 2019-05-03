% Main script to run other scripts for theta model v5

clear; close all; clc

tic

% A is connectivity matrix (1-ER, 2-small world, 3-scale free, 4-hierarchy)
A = 3;
<<<<<<< HEAD:v7/ablationstudy/main.m
ablate = 1;
=======
ablate = 0;
>>>>>>> 3bb8e1c918cc3939724943d01d1e4ff33bf16595:v7/ablationstudy_networkstudy/main.m

synctheta_v7_edit(6e4,1,1,0,A,ablate);

toc
