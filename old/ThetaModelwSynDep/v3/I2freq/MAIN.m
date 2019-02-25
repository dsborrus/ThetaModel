% Script to map the I parameter to the idea of a neuron's intrinsic
% frequency
clear; close all; clc; tic;
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

i_min = -.01;
i_max = 2;
resolution = 10;


iu = linspace(i_min,i_max,resolution);

f_array = zeros(resolution,1);

for k = 1:resolution
    f_array(k) = synctheta_I2F(iu(k));
    toc
    disp(['Finished ' mat2str(k) ' out of ' mat2str(resolution)]);
end

%%
figure
plot(iu,f_array); hold on;
xlabel('I')
ylabel('f (spikes/second)')
title('Intrinsic frequency vs. I')
    
