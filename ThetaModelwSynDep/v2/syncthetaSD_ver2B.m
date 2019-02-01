% Sync Theta Model with synaptic depression - Version 1 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% Original ode was simple theta model This now has synaptic depresseion
%
% Version 2 - no longer going to store massive theta and y array
% just a few key neurons and spike average
% Changed theta(j,:) and theta(j+1,:) to thetaNOW and thetaNEXT
% respectively. To try and save time. Didn't help much...
% Also changed the way spikes is calculated. Before ss would compute how
% many spikes all the neurons RECEIEVED at that time step. However, I think
% it's better to have ss calculate all the spikes the neurons SENT in that
% time step. This will always be less than n, unlike before.
%
% Version 1 - including synaptic depression
% 
% System of linked theta neurons, maybe with two populations of different
% intrinsic frequency
% Also there is a synaptic depression
%
% Function is
% d theta/ d time = 1 - cos(theta) + [1 + cos(theta)] * I 
% d y / d t = (y - 1)/tau where on every burst, y = y-ydrop
% hint: positive I is instrinicly rhymic

clc
clear
close all
tic

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

% % % User Params % % %

n1 = 100;   % number of neurons in the first population
n2 = 0;     % number of neurons in the second population
dt = 0.01;   % time step
tmax = 5e2;   % maximum time for simulation
iu1 = 0;  % mean I parameter for first population
isig1 = 0.1;  % std of I parameter for first population
iu2 = 0.0001;  % mean I parameter for second population# 
isig2=0.00001;    % std of I parameter for second population
prob = 0.85; % E-R graph, prob is prob of connection.
D=3;      % Strength of networkness
tauavg=1;   % Relaxation of network excitement

ydrop = .2; % How much of an effect firing has on synaptic depression
            % (should be between 0 and 1)!!!
tauy  =  10; % Char time for return to ss for y (synap depress)         

% Which neurons to save
impNeu = [1, 3, 5];
           

% % % Script Settings % % %

DoDBPlot = 1; % Should we make the graph dan made?

% % % % Script % % % %

n=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
tnum = tmax/dt; % the nummber of time steps
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps
phase = pmin:0.01:pmax; 

% set intialize values for theta
thetaNOW = pmin + (pmax-pmin).*rand(1,n);
% set initial values for y
yNOW = (randn(1,n)*.001)+.5;
% initialize the ODE future point output array
thetaNEXT = NaN(n,1);
yNEXT = NaN(n,1);
% initialize ODE output for key neurons
thetaIMP = zeros(tnum,length(impNeu));
yIMP     = zeros(tnum,length(impNeu));
% set first value
thetaIMP(1,:) = thetaNOW(impNeu);
yIMP(1,:)     = yNOW(impNeu);

% initialize the spike array
spikes = NaN*ones(tnum,1); spikes(1)=0;

%rplot = NaN*ones(tnum,1); rplot(1)=0;   % initialize the raster plot?

% initialize I vector recall, it should be length of n1+n2
I = [ iu1+isig1*randn(1,n1) iu2+isig1*randn(1,n2) ] 

% strength of connections (total excitability of network 
% divided by the number of connections)
delta = D/n; 

% creates connectivity matrix, damn that's savy
A=(rand(n,n)<prob);
Atrans = A'; % time saving thing

% makes sure no neuron is connected to itself. Again, quite savy
for m=1:n, A(m,m)=0; end

% note: gif making feature has been removed

% Begin simulation loop
for j = 1:tnum-1
    
    % Calculate ODEs next step (Euler's method)
    thetaNEXT = thetaNOW + dt * thetaODE(thetaNOW,I);
    yNEXT     = yNOW     + dt * yODE(yNOW,tauy);
    
    
    e=1; ss=0;
    while any(e)
    % check if any neuron is above pi ( should be pmax )
    e = (thetaNEXT>pi); a=find(e);
    % reset that neuron back 2pi
    thetaNEXT(a)=thetaNEXT(a)-2*pi;
    % augment the synpatic depression term for next y
    yNEXT(a) = yNEXT(a) - ydrop;
    % calculate the number of pulses the downstream neuron receives
    %s_old=e*A; % but now get deeper....
    % calculate the strength + number of the pulses (1 is max strength), though 
    % HERE'S A QUESTION THAT NEEDS TO BE RESOLVED
    % SHOULD I MULTIPLY BY y(j) OR y(j+1)??? I THINK y(j). Because it is
    % the value of synaptic depression before the current spike that
    % matters?
    % anyways, first I multiple e element wise with the synaptic strength
    % of each neuron, before multiplying by A
    e_power = e.*yNOW;
    
    s = (Atrans*e_power')';
    
    % Calculate and add the pulse
    thetaNEXT = thetaNEXT+delta*s;
    
    % Sum up number of spikes in this step
    %ss = ss+sum(s_old);
    ss = ss+sum(length(a));
    
    % little check for sanity %
    if any(delta*s>2*pi)
        error('I think this means we pushed a neuron up 2pi, it missed a phase')
    end
    
    end

    spikes(j+1) = spikes(j)+dt*(ss/n-spikes(j)/tauavg);
    %rplot(j+1) = 1/n*sum(exp(i*theta(j+1,:)));
    
    % Also save a few key neurons
    thetaIMP(j+1,:) = thetaNEXT(impNeu);
    yIMP(j+1,:)     = yNEXT(impNeu);
    
    % Now move Next back to Now
    thetaNOW = thetaNEXT;
    yNOW = yNEXT;
    
end

if DoDBPlot
    DBPlot
end

toc

function dtheta = thetaODE(theta,I)
    dtheta = 1-cos(theta) + (1+cos(theta)).*I;
end

function dy = yODE(y,tauy)
    dy = (1 - y)/tauy;
end


