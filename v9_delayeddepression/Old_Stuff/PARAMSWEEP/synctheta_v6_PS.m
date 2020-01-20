function synctheta_v6_PS(tmax,params)
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
%
% Version 5 - updating pulses, no longe rpulse coupled, but upstream
% neurons create synaptic current in downstream neuron
% Notes on v5: it is trash. Could never get it to oscillate correctly :(
% v6 will drop the continous synaptic current and use discontinous
% like synaptic depression is
% I think the problem was in the histerises loop? See FindingBistability
% folder (that's actually the only folder that got changed along with main
% directoty.
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

n1 = params.n1;
n2 = 0;     % number of neurons in the second population
dt = 0.05;   % time step
%tmax = 1e4;     % maximum time of simulation
iu1 = 0;  % mean I parameter for first population
isig1 = params.isig;  % std of I parameter for first population
iu2 = 0;  % mean I parameter for second population#
isig2 = 1e-6;    % std of I parameter for second population
prob = params.prob; % E-R graph, prob is prob of connection.
D = params.D;      % Strength of networkness
tauavg = 1e2;   % Relaxation of network excitement

ydrop = params.ydrop; % How much of an affect firing has on synaptic depression
% (should be between 0 and 1)!!!
tauy  =  5e3; % Char time for return to ss for y (synap depress)

sigain = 1.0;
tausi = params.tausig;

istate=4;

% % % Script Settings % % %

DoDBPlot = 1; % Should we make the graph dan made?

% % % % Script Stuff % % (No need to edit below this line if casual user)

n=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
tnum = tmax/dt; % the nummber of time steps
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps

% initialize the ODE output arrays(theta array and y array)
theta = zeros(tnum,n);
y     = zeros(tnum,n);
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

disp(['Intrinsic freq that are positive = ' mat2str(length(find(I>0))*100/length(I)) '%'])

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
        y(1,:)     = (randn(1,n)*.01) + 0.1;
        %sij(1,:)
    case 2 %high state
        theta(1,:) = pmax.*rand(1,n);
        y(1,:)     = (randn(1,n)*.01) + 0.97;
        si(1,:)   = (randn(1,n)*.01) + 0.5;
    case 3 %med state
        theta(1,:) = pmin + (pmax-pmin).*rand(1,n);
        y(1,:) = (randn(1,n)*.07)+.8;
        si(1,:)   = (randn(1,n)*.2) + 0.5;
    case 4
        theta(1,:) = pmin + (pmax-pmin).*rand(1,n);
        y(1,:) = ones(n,1);
        si(1,:)   = (randn(1,n)*.2) + 0.5;
end

wb = waitbar(0,'Simulating Simulation');

% Begin simulation loop
for j = 1:tnum-1
    
    if mod(j,2000) == 0
        waitbar(j/tnum,wb)
    end
    
    % calculate synaptic strengths
    % Sum of incoming current is = 
    % instrinsic current (I) +
    % strength of connections (delta) *
    % synaptic depression (y) *
    % synaptic current generated (sij) matrix mulitplied by
    % connectivity matrix (A)
    Isummed = I + (delta*y(j,:).*si(j,:)) * A; 
    %Isummed2 = I + ( sum(A(e,:).*(delta*y(j,e).*sij(j,:))',1));
    %Isummed = I + ( sum(A(e,:).*y(j,e)',1).*sij(j,:));
    Ihistory(j) = Isummed(rr);
    sihistory(j) = si(j,:) * A(:,rr);

    % Calculate ODEs next step (Euler's method)
    theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),Isummed);
    y(j+1,:)     = y(j,:)     + dt * yODE(y(j,:),tauy);
    si(j+1,:)   = si(j,:)   + dt * siODE(si(j,:),tausi);
    
    e=1; ss=0;
    while any(e)
        % check if any neuron is above pi ( should be pmax )
        e = (theta(j+1,:)>pi); a=find(e);
        % reset that neuron back 2pi
        theta(j+1,a)=theta(j+1,a)-2*pi;
        % augment the synpatic depression term for next y
        y(j+1,a) = y(j+1,a) - ydrop*y(j+1,a);
        si(j+1,a) = si(j+1,a) + sigain*(1-si(j+1,a));
        
        % count spikes
        ss = ss+length(a);
        
        % catalouge for raster plot
        raster(j,a) = 1;
        
    end
    
    spikes(j+1) = spikes(j)+ss/n+dt*(-spikes(j)/tauavg);
    
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
    vin(11) = sigain;
    vin(12) = tausi;
    
    DBPlot_v2(dt,tmax,t,n,y,si,spikes,raster,rr,Ihistory,sihistory,vin,tauavg,istate)
end

close(wb)


% ODE functions

    function dtheta = thetaODE(theta,I)
        dtheta = 1-cos(theta) + (1+cos(theta)).*I;
    end

    function dy = yODE(y,tauy)
        yss = 1;
        dy = (yss - y)/tauy;
    end

    function dsi = siODE(si,tausi)
        siss = 0;
        dsi = (siss-si)/tausi;
    end

end