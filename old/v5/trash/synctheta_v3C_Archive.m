% Sync Theta Model with synaptic depression - Version 3 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% Original ode was simple theta model This now has synaptic depresseion
%
% % % % Update History % % % %
%
% Miniupdate C - Increase speed of simulation by removing the large matrix
% multiplication for connectivity matrix. Now picks and chooses only
% neurons it needs!
%
% Version 3-  Update name: Restore.
% Now, changing the theta pulse. Before it was D/n.
% But we must also accoount for the connectivity of the network
% D/((N-1)P).
% Also graphics are getting an update. Need raster plot (DB plot)
%
%
% Version 2 - Tried to speed up code, it failed. Restored back to version 1
%
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

% % % Begin % % %

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
tmax = 1e3;     % maximum time of simulation
iu1 = -0.01;  % mean I parameter for first population
isig1 = 0.00;  % std of I parameter for first population
iu2 = 0.001;  % mean I parameter for second population# 
isig2 = 0.0001;    % std of I parameter for second population
prob = 0.85; % E-R graph, prob is prob of connection.
D = 1;      % Strength of networkness
tauavg=1;   % Relaxation of network excitement

ydrop = .2; % How much of an effect firing has on synaptic depression
            % (should be between 0 and 1)!!!
tauy  =  6; % Char time for return to ss for y (synap depress)           
           

% % % Script Settings % % %

nplot=5;   % How long to wait before plotting gif
DoDBPlot = 0; % Should we make the graph dan made?

% % % % Script Stuff % % (No need to edit below this line if casual user)

n=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
tnum = tmax/dt; % the nummber of time steps
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps
phase = pmin:0.01:pmax; 

% initialize the ODE output arrays(theta array and y array)
theta = zeros(tnum,n);
y = zeros(tnum,n);

spikes = NaN*ones(tnum,1); spikes(1)=0; % initialize the spike array

%rplot = NaN*ones(tnum,1); rplot(1)=0;   % initialize the raster plot?

% initialize I vector recall, it should be length of n1+n2
I = [ iu1+isig1*randn(1,n1) iu2+isig1*randn(1,n2) ]; 

% strength of connections (total excitability of network 
% divided by the number of connections)
delta = D/n; 

% creates connectivity matrix, damn that's savy
A=(rand(n,n)<prob);

% makes sure no neuron is connected to itself. Again, quite savy
for m=1:n, A(m,m)=0; end

% set intialize values for theta
theta(1,:) = pmin + (pmax-pmin).*rand(1,n);
% set initial values for y
y(1,:) = (randn(1,n)*.001)+.5;

% note: gif making feature has been removed

% Begin simulation loop
for j = 1:tnum-1
    
    % Calculate ODEs next step (Euler's method)
    theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),I);
    y(j+1,:)     = y(j,:)     + dt * yODE(y(j,:),tauy);
    
    
    e=1; ss=0;
    while any(e)
    % check if any neuron is above pi ( should be pmax )
    e = (theta(j+1,:)>pi); a=find(e);
    % reset that neuron back 2pi
    theta(j+1,a)=theta(j+1,a)-2*pi;
    % augment the synpatic depression term for next y
    y(j+1,a) = y(j+1,a) - ydrop;
    % calculate the number of pulses the downstream neuron receives
    %s_old=e*A; % but now get deeper....
    % calculate the strength + number of the pulses (1 is max strength), though 
    % HERE'S A QUESTION THAT NEEDS TO BE RESOLVED
    % SHOULD I MULTIPLY BY y(j) OR y(j+1)??? I THINK y(j). Because it is
    % the value of synaptic depression before the current spike that
    % matters?
    % anyways, first I multiple e element wise with the synaptic strength
    % of each neuron, before multiplying by A
    %e_power = e.*y(j,:);
    
    %s = e_power*A;
    %Apower1 = A.*y(j,:)';
    s = sum(A(e,:).*y(j,e)');
    %s1 = e_power*A;
    %s2 = e*Apower;
    %s1 = sum(Apower1(e,:),1);
    
    % Calculate and add the pulse
    theta(j+1,:) = theta(j+1,:)+delta*s;
    
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
