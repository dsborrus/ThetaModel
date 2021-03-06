function [sps] = synctheta_v6_PS(D,params)
%
% code for parameter sweep. Taking out synaptice depression.
%
%
% Sync Theta Model with synaptic depression - Version 5 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% Original ode was simple theta model This now has synaptic depresseion
%
% tmax is maximum time to run simulation
% istate is initial condition (1-low, 2-high, 3-random spread)
%
% % % % Update History % % % %
%
% Version 6 - removing continous synaptic current, replacing it with
% discontinous jumps, super similar to synpatic depression
%
% Version 5 - updating pulses, no longe rpulse coupled, but upstream
% neurons create synaptic current in downstream neuron
%
% Version 4 - allowing for multiple spikes in one time "frame"
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

% unit for time is ms

% % % User Params % % %

% duration of bump (if we are doing one)
% dur = 30;

n1 = params.n1;   % number of neurons in the first population
n2 = params.n2;     % number of neurons in the second population
dt = params.dt;   % time step
tmax = params.tmax; % maximum time of simulation
iu1 = params.iu1;  % mean I parameter for first population
isig1 = params.isig1;  % std of I parameter for first population
iu2 = params.iu2;  % mean I parameter for second population#
isig2 = params.isig2;    % std of I parameter for second population
prob = params.prob; % E-R graph, prob is prob of connection.
%D = D;      % Strength of networkness
tauavg=params.tauavg;   % Relaxation of network excitement

sigain = params.sigain;
tausi = params.tausi;

% 
% ydrop = 0; % How much of an affect firing has on synaptic depression
% % (should be between 0 and 1)!!!
% tauy  =  5e3; % Char time for return to ss for y (synap depress)

silence = params.silence;
bumpit = params.bumpit;
istate = params.istate;

% % % Script Settings % % %

DoDBPlot = 0; % Should we make the graph dan made?

% % % % Script Stuff % % (No need to edit below this line if casual user)

n=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
tnum = tmax/dt; % the nummber of time steps
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps

% initialize the ODE output arrays(theta array and y array)
theta = zeros(tnum,n);
y     = ones(tnum,n);
si   = zeros(tnum,n);

% initialize auxiallry outputs
spikes = NaN(tnum,1); spikes(1)=0; % initialize the spike array
raster = NaN(tnum,n); % initialize the raster plot
Ihistory = NaN(tnum,1);
sihistory = NaN(tnum,1);
% random neuron to track
rr = randi(n1);

% initialize I vector recall, it should be length of n1+n2
I = [ iu1+isig1*randn(1,n1) iu2+isig2*randn(1,n2) ];

% strength of connections (total excitability of network
% divided by the number of neurons (minus 1) and divided by probability
% of connections)
delta = D/((n-1)*prob);

% creates connectivity matrix, damn that's savy
A=(rand(n,n)<prob);

% makes sure no neuron is connected to itself. Again, quite savy
for m=1:n, A(m,m)=0; end

% Set initial values for theta and y
switch istate
    case 1 %low state
        theta(1,:) = pmin+randn(1,n) * .01;
        %y(1,:)     = (randn(1,n)*.01) + 0.1;
        %sij(1,:)
    case 2 %high state
        theta(1,:) = pmax.*rand(1,n);
        %y(1,:)     = (randn(1,n)*.01) + 0.97;
        si(1,:)   = (randn(1,n)*.01) + 0.5;
    case 3 %med state
        theta(1,:) = pmin + (pmax-pmin).*rand(1,n);
        %y(1,:) = (randn(1,n)*.07)+.8;
        si(1,:)   = (randn(1,n)*.2) + 0.5;
end

% initialize spike counter
% spikes per frame
spf = zeros(tnum,1);

% Begin simulation loop
for j = 1:tnum-1
    
    % enforce silent period to allow equillibrium
	if silence
       if t(j) < 1000
          si(j,:) = 0; 
       end
    end
    
    % calculate synaptic strengths
    % Sum of incoming current is = 
    % instrinsic current (I) +
    % strength of connections (delta) *
    % synaptic depression (y) *
    % synaptic current generated (sij) matrix mulitplied by
    % connectivity matrix (A)
    %Isummed = I + (delta*y(j,:).*si(j,:)) * A; 
    Isummed = I + (delta*si(j,:)) * A; 
    
    % bumping system to boot start it
    if bumpit
        if t(j) > 500 && t(j) < 1000
                Isummed = Isummed + 0.5;
        elseif t(j) > 1500 && t(j) < 2000
                Isummed = Isummed + 0.5;
        elseif t(j) > 2500 && t(j) < 3000
                Isummed = Isummed + 0.5;
        end
    end
    
    
    % record keeping
    Ihistory(j) = Isummed(rr);
    sihistory(j) = si(j,:) * A(:,rr);

    % Calculate ODEs next step (Euler's method)
    theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),Isummed);
    %y(j+1,:)     = y(j,:)     + dt * yODE(y(j,:),tauy);
    si(j+1,:)   = si(j,:)   + dt * siODE(si(j,:),tausi);
    
    e=1; ss=0;
    while any(e)
        % check if any neuron is above pi ( should be pmax )
        e = (theta(j+1,:)>pi); a=find(e);
        % reset that neuron back 2pi
        theta(j+1,a)=theta(j+1,a)-2*pi;
        % augment the synpatic depression term for next y
        %y(j+1,a) = y(j+1,a) - ydrop*y(j+1,a);
        si(j+1,a) = si(j+1,a) + sigain*(1-si(j+1,a));
        
        % count spikes
        ss = ss+length(a);
        
        % catalouge for raster plot
        raster(j,a) = 1;
        
    end
    
    spikes(j+1) = spikes(j)+dt*(ss/n-spikes(j)/tauavg);
    
    spf(j) = ss;
end

% plotting

if DoDBPlot
    
    vin(1) = D;
    vin(2) = isig1;
    vin(3) = iu1;
    vin(4) = isig2;
    vin(5) = iu2;
    vin(6) = n;
    vin(7) = prob;
    vin(8) = tauavg;
    vin(9) = tauy;
    vin(10) = ydrop;
    
    DBPlot_v2(dt,tmax,t,n,y,si,spikes,raster,rr,Ihistory,sihistory,vin)
end

%% Analysis %%

% time for analysis window
tstart_an = tmax/2;
tend_an = tmax-1;

tw_an = tstart_an/dt:tend_an/dt;

spf_insidewindow = spf(tw_an);

sps = sum(spf_insidewindow)/(length(spf_insidewindow)*dt/1000);

%% ODE functions

    function dtheta = thetaODE(theta,I)
        dtheta = 1-cos(theta) + (1+cos(theta)).*I;
    end

    function dy = yODE(y,tauy)
        dy = (1 - y)/tauy;
    end

    function dsi = siODE(si,tausi)
        siss = 0;
        dsi = (siss-si)/tausi;
    end

end