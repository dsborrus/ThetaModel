% Main script to run run param sweep

clear; close all; clc

tic

% % % User Params % % %
n1_array     = [100;   200;    500];
prob_array   = [1.0];
D_array      = [0.001; 0.0008; 0.0005; 0.0004; 0.0001];
ydrop_array  = [0.5;   0.3;    0.1;    0.05];
tausig_array = [1;     5;      10;     20];
isig_array   = [1e-6;  1e-7;   1e-8];

num_sims = length(n1_array)*length(prob_array)*length(D_array)*length(...
    ydrop_array)*length(tausig_array)*length(isig_array)

waitforbuttonpress

for i_n1 = 1:length(n1_array)
for i_pr = 1:length(prob_array)
for i_D  = 1:length(D_array)
for i_yd = 1:length(ydrop_array)
for i_ta = 1:length(tausig_array)
for i_is = 1:length(isig_array)    
    params.n1 = n1_array(i_n1);
    params.prob = prob_array(i_pr); 
    params.D = D_array(i_D);     
    params.ydrop = ydrop_array(i_yd);
    params.tausig = tausig_array(i_ta);
    params.isig = isig_array(i_is);
    synctheta_v6_PS(3e4,params)
end
end
end
end
end
end
    
istate=4;

synctheta_v6_PS(5e4,4)


toc

