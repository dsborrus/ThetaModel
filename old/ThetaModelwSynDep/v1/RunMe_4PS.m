% Main Script for parameter sweep
% written by DB 1/24/19
clc
clear all
close all
tic
 %#ok<*NBRAK>

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

n1 = [100 120];   % number of neurons in the first population
n2 = [0];     % number of neurons in the second population
iu1 = [0.001 0.01];  % mean I parameter for first population
isig1 = [0.001 0.01];  % std of I parameter for first population
iu2 = [2.0];  % mean I parameter for second population
isig2 = [0.0];    % std of I parameter for second population
prob = [0.9 0.85]; % E-R graph, prob is prob of connection.
D = [0.1 0.05];      % Strength of networkness
tauavg = [0.5 1 2 5];   % Relaxation of network excitement

iter = 0;

MaxIter = length(n1)*length(n2)*length(iu1)*length(isig1)*length(iu2)*...
    length(isig2)*length(prob)*length(D)*length(tauavg);

for in1 = 1:length(n1)
for in2 = 1:length(n2)
for iiu1 = 1:length(iu1)
for iisig1 = 1:length(isig1)
for iiu2 = 1:length(iu2)
for iisig2 = 1:length(isig2)
for iprob = 1:length(prob)
for iD = 1:length(D)    
for itauavg = 1:length(tauavg)
    
    %[gradeMa(in1,in2,iiu1,iisig1,iiu2,iisig2,iprob,iD,itauavg),...
    %    gradeMu(in1,in2,iiu1,iisig1,iiu2,iisig2,iprob,iD,itauavg)] = ...
    %    synctheta_ver6PS(n1(in1),n2(in2),iu1(iiu1),isig1(iisig1),iu2(iiu2),...
    %    isig2(iisig2),prob(iprob),D(iD),tauavg(itauavg));
    
    [gradeMa,gradeMu] = ...
        synctheta_ver6PS(n1(in1),n2(in2),iu1(iiu1),isig1(iisig1),iu2(iiu2),...
        isig2(iisig2),prob(iprob),D(iD),tauavg(itauavg));
    iter = iter + 1;
    gradebook.page(iter).params = [n1(in1),n2(in2),iu1(iiu1),isig1(iisig1)...
        ,iu2(iiu2),isig2(iisig2),prob(iprob),D(iD),tauavg(itauavg)];
    gradebook.page(iter).GMax = gradeMa;
    gradebook.page(iter).GMu  = gradeMu;
    
    
    %disp(['tested ' mat2str([in1,in2,iiu1,iisig1,iiu2,iisig2,iprob,iD,itauavg])])
    disp(['Completed ' mat2str(iter) ' of ' mat2str(MaxIter) ' combinations'])
    
end
end
end
end
end
end
end
end
end

[~,i] = maxk([gradebook.page(:).GMax],3)
gradebook.page(i(1)).params


toc
