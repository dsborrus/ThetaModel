function synctheta_v7(tmax,istate,p1,p2)
%
% Sync Theta Model with synaptic depression - Version 5 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% Original ode was simple theta model This now has synaptic depresseion
%
% tmax is maximum time to run simulation
% istate is initial condition (1-low, 2-high, 3-random spread, 4 - previous state)
%
% p1 is the main plot (1 for yes. 0 for no)
% p2 is the gif-animated plot (1 for yes. 0 for no) recommend most users
% for p2 = 0
%
%
% % % % Update History % % % %
%
% Version 7 - delaying synaptic depression
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

n1 = 100;   % number of neurons in the first population
n2 = 0;     % number of neurons in the second population
dt = 0.2;   % time step
%tmax = 1e4;     % maximum time of simulation
iu1 = -.0009;  % mean I parameter for first population
isig1 = 1e-7;  % std of I parameter for first population
iu2 = 0;  % mean I parameter for second population#
isig2 = 0;    % std of I parameter for second population
prob = .5; % E-R graph, prob is prob of connection.
D = 0.03;      % Strength of networkness
tauavg = 1e2;   % Relaxation of network excitement

tautheta = 1;

mgain = .2; % How much of an affect firing has on synaptic depression
taum  =  300; % Char time for return to ss for m (synap depress)

nrise = 0.011;
taun  = 1300;

sigain = 1;
tausi = 15;

noisesigma=.0090;

% % % Script Settings % % %

DoDBPlot = p1; % Should we make the graph dan made?
DoPDPlot = p2; % Should we make the gif of the Phase Diagrams of all neurons?

% % % % Script Stuff % % (No need to edit below this line if casual user)

N=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
tnum = tmax/dt; % the nummber of time steps
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps

% initialize the ODE output arrays(theta array and y array)
theta = zeros(tnum,N);
y     = zeros(tnum,N);
si    = zeros(tnum,N);
m     = zeros(tnum,N);
n     = zeros(tnum,N);

% initialize auxiallry outputs
spikes = NaN(tnum,1); spikes(1)=0; % initialize the spike array
raster = NaN(tnum,N); % initialize the raster plot
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
delta = D/((N-1)*prob);

% creates connectivity matrix, damn that's savy
A=(rand(N,N)<prob);

% makes sure no neuron is connected to itself. Again, quite savy
for i=1:N, A(i,i)=0; end

% Set initial values for theta and y
switch istate
    case 1 %low state
        theta(1,:) = pmin+randn(1,N) * .01;
        y(1,:)     = zeros(1,N) + 0.1;
        si(1,:)    = zeros(1,N) + 0.1;
        m(1,:)     = zeros(1,N);
        n(1,:)     = zeros(1,N);
    case 2 %high state
        theta(1,:) = pmax.*rand(1,N);
        y(1,:)     = ones(1,N);
        si(1,:)    = zeros(1,N) + 0.5;
        m(1,:)     = zeros(1,N)+1;
        n(1,:)     = zeros(1,N)+1;        
    case 3 %med state
        theta(1,:) = pmin + (pmax-pmin).*rand(1,N);
        y(1,:)     = randn(1,N)*.1 + 0.8;
        si(1,:)    = randn(1,N)*.2 + 0.5;
        m(1,:)     = randn(1,N)*.2 + 0.5;
        n(1,:)     = zeros(1,N)*.2 + 0.5;  
    case 4 %last state of previous simulation
        load lastconditions lc
        if lc.N == N
        theta(1,:) = lc.theta;
        y(1,:)     = lc.y;
        si(1,:)    = lc.si;
        m(1,:)     = lc.m;
        n(1,:)     = lc.n; 
        else
        error('cant use last conditions, number of neurons changed')
        end
end

%% graphics stuff please skip
wb = waitbar(0,'Simulating Simulation');
gifresolution = 10; %ms
if DoPDPlot
    gf = figure('position',[10 100 2200 1100]);
    subplot(6,6,[1 21]); hold on;
    axis tight manual; hold on
    ylim([-.01 .01]);
    xlim([-0.17 0.17]);
    gffilename = 'testAnimated.gif';
    plot([pmin pmax],zeros(2,1),'k--','linewidth',0.5); hold on;
    P = []; P2 = []; P3 = [];
    title('Mean I of network')
    
    subplot(6,6,[31 36]);
    hold on

    str1 = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
        '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
        '. N = ' mat2str(N) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
        '. mgain = ' mat2str(mgain) '. taum = ' mat2str(taum) ...
        '. sigain = ' mat2str(sigain) '. tausi = ' mat2str(tausi)];
    annotation('textbox',[.1 .08 .1 .1],'String',str1,'FitBoxToText','on');

    str2 = ['nrise = ' mat2str(nrise) '. taun = ' mat2str(taun) ...
           '. tautheta = ' mat2str(tautheta) '. noisesigma = ' mat2str(noisesigma)...
           '. istate = ' mat2str(istate)];
    annotation('textbox',[.1 .03 .1 .1],'String',str2,'FitBoxToText','on');


    axis off
end

vin(1) = D;
vin(2) = isig1;
vin(3) = iu1;
vin(4) = isig2;
vin(5) = iu2;
vin(6) = N;
vin(7) = prob;
vin(8) = tauavg;
vin(9) = mgain;
vin(10) = taum;
vin(11) = nrise;
vin(12) = taun;
vin(13) = sigain;
vin(14) = tausi;
vin(15) = tautheta;
vin(16) = noisesigma;

%% Begin simulation loop, come back now
for j = 1:tnum-1
    
    if mod(j,2000) == 0
        waitbar(j/tnum,wb)
    end
    
    % calculate synaptic strengths
    Isummed = I + ( (delta*y(j,:).*si(j,:)) * A); 
    
    % Record this
    Ihistory(j) = Isummed(rr);
    sihistory(j) = si(j,:) * A(:,rr);

    % Calculate ODEs next step (Euler's method)
    theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),Isummed,tautheta)...
                   + noise(dt,N,noisesigma)';
    m(j+1,:) = m(j,:) + dt * mODE(m(j,:),taum);
    n(j+1,:) = n(j,:) + dt * nODE(n(j,:),m(j,:),nrise,taun);
    si(j+1,:)   = si(j,:)   + dt * siODE(si(j,:),tausi);
    
    y(j+1,:) = 1 - n(j+1,:);
    
    e=1; ss=0;
    while any(e)
        % check if any neuron is above pi ( should be pmax )
        e = (theta(j+1,:)>pi); a=find(e);
        % reset that neuron back 2pi
        theta(j+1,a)=theta(j+1,a)-2*pi;
        % augment the synpatic depression term & conductance for next j
        m(j+1,a) = m(j+1,a) + mgain * (1-m(j+1,a));
        si(j+1,a) = si(j+1,a) + sigain*(1-si(j+1,a));
        
        % count spikes
        ss = ss+length(a);
        
        % catalouge for raster plot
        raster(j,a) = 1;
        
    end
    
    spikes(j+1) = spikes(j)+ss/N+dt*(-spikes(j)/tauavg);
    
    
    %%% giffing %%%%
    if mod(j-1,20*gifresolution) == 0 && DoPDPlot
        if j > 20*3000 % over 2 seconds
            [P,P2,P3] = updategif(theta(j,:),pmin,pmax,Isummed,gf,gffilename,j,P,P2,P3,raster,spikes,t,vin,istate);
        end
    end
    
end

% plotting

if DoDBPlot
    DBPlot_v2(dt,tmax,t,y,m,n,si,spikes,raster,rr,Ihistory,sihistory,vin,istate)
end

close(wb)

% save last conditions
lc.theta = theta(end,:);
lc.y = y(end,:);
lc.m = m(end,:);
lc.n = n(end,:);
lc.si = si(end,:);
lc.N = N;

save('lastconditions.mat','lc')


% save (send) output
%thetas = theta;


% ODE functions

    function dtheta = thetaODE(theta,I,tautheta)
        dtheta = (1/tautheta)* (1-cos(theta)) + (1+cos(theta)).*I;
    end

    function dm = mODE(m,taum)
        dm = -m/taum;
    end

    function dn = nODE(n,m,nrise,taun)
        dn = nrise*m.*(1-n) - n/taun;
    end

    function dsi = siODE(si,tausi)
        siss = 0;
        dsi = (siss-si)/tausi;
    end

    function noiseout = noise(deltat,NumNeurons,sigma)
    
        noiseout = sigma * sqrt(deltat) * randn(NumNeurons,1);
        
    end

end