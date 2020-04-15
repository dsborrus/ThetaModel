function outputs = simulate_v1(Parameters,Options)
%
% Sync Theta Model with synaptic depression - Version 8 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% ConMat matrix done by Rachel Smith
%
% tmax is maximum time to run simulation
% istate is initial condition (1-low, 2-high, 3-random spread, 4 - previous state)
%
%
% This is the function, run the a main* script to call this function
%
% Update history moved to README
%
% unit for time is ms

%% Load Parameters
disp('Loading parameters')

% System params %
tmax     = Parameters.tmax;    % maximum time of simulation
dt       = Parameters.dt;      % time step

% Neurons params %

iu1      = Parameters.iu1;     % mean I parameter for first population
isig1    = Parameters.isig1;   % std of I parameter for first population
iu2      = Parameters.iu2;     % mean I parameter for second population#
isig2    = Parameters.isig2;   % std of I parameter for second population
tautheta = Parameters.tautheta;% relaxation of neuron's theta
    % synaptic depression 1
mgain    = Parameters.mgain;   % How much of an affect firing has on synaptic depression
taum     = Parameters.taum;    % Char time for return to ss for m (synap depress)
    % synaptic depression 2
nrise    = Parameters.nrise;   % Rate of rise for synaptic depression based on n
taun     = Parameters.taun;    % Char time for return to ss for n (synap depress)
    % snynaptic conductance
sigain   = Parameters.sigain;  % How much of a gain firing has on synaptic conductance
tausi    = Parameters.tausi;   % Char time for return to ss for n (synap depress)
    % noise
noisesig = Parameters.noisesig;% Variance of noise

% Network params %

n1       = Parameters.n1;      % number of neurons in the first population
n2       = Parameters.n2;      % number of neurons in the second population
D        = Parameters.D;       % Strength of networkness
tauavg   = Parameters.tauavg;  % Relaxation of network excitement
istate   = Parameters.istate;  % initial state to use
conmat   = Options.conmat;     % network architecture to use
    % case 1 - ER
prob     = Parameters.er_prob;    % prob of connection
    % case 2 - small world
sw_M     = Parameters.sw_M;    % number of Ns on each side
sw_p     = Parameters.sw_p;    % probability of "short cut" 
    % case 3 - scale-free
sf_mo   = Parameters.sf_mo;    % size of seed
sf_m    = Parameters.sf_m;     % average degree (use mo=m or m<mo)

% ablation params %
Doablate = Options.Doablate;
t2ablat = Parameters.t2ablat;  % When to ablate neurons (every XXX seconds)
N2k_w1p = Parameters.N2k_w1p;  % number of neurons to kill with one pulse

% Script params %
doAplot = Options.doAplot;
doplot1 = Options.doplot1;
doplot2 = Options.doplot2;
dogifplot = Options.dogifplot;

%% Script Stuff (No need to edit below this line if casual user)
disp('Pre-sim stuff...')

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
if Options.trackstatevariablesmeans
    meany = zeros(tnum,1);
    meansi = zeros(tnum,1);
end


% random neuron to track
rr = randi(n1);

% initialize I vector recall, it should be length of n1+n2
I = [ iu1+isig1*randn(1,n1) iu2+isig2*randn(1,n2) ];

disp([ mat2str(sum(I>0)) ' out of ' mat2str(length(I)) ' neurons are tonic'])

%% Set initial values for theta and y
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

if Options.trackstatevariablesmeans
    meany(1) = mean(y(1,:));
    meansi(1) = mean(si(1,:));
end

%%  make connectivity matrix
disp('Making connectivity matrix')
A = MakeNetwork(Parameters,Options);

%% more pre-sim stuff
% strength of connections (total excitability of network
% divided by the number of neurons (minus 1) and divided by probability
% of connections)
%delta = D/N/sum(sum(A));
meanInDegree = mean(sum(A));
delta = D/meanInDegree;

%initializing ablation scheme (if doing it)
killlist = randperm(100);
k_indx = 0;
getnoisy = true(N,1);

%% gif graphics stuff (and waitbar) 
wb = waitbar(0,'Simulating Simulation');
gifresolution = 10; %ms
if dogifplot
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
           '. tautheta = ' mat2str(tautheta) '. noisesigma = ' mat2str(noisesig)...
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
vin(16) = noisesig;

%% Begin simulation loop
disp('Beginning simulation')
for j = 1:tnum-1
    
    if mod(j,2000) == 0
        waitbar(j/tnum,wb)
    end
    
    %laser ablation
    if mod(j, t2ablat * 1000/dt) == 0 && Doablate
        for lala = 1:N2k_w1p
            k_indx = k_indx + 1;
            deadN = killlist(k_indx);
            A(deadN, :) = 0;
            A(:, deadN) = 0;
            % update get noisy, dead neurons dont get noise
            getnoisy(deadN) = false;
        end
    end
    
    % calculate synaptic strengths
    Isummed = I + ( (delta*y(j,:).*si(j,:)) * A); 
    
    % Record this
    Ihistory(j) = Isummed(rr);
    sihistory(j) = si(j,:) * A(:,rr);

    % Calculate ODEs next step (Euler's method)
    theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),Isummed,tautheta)...
                   + noise(dt,N,noisesig,getnoisy)';
    m(j+1,:) = m(j,:) + dt * mODE(m(j,:),taum); 
    n(j+1,:) = n(j,:) + dt * nODE(n(j,:),m(j,:),nrise,taun);
    si(j+1,:)   = si(j,:)   + dt * siODE(si(j,:),tausi);
    
%     % turning off synaptic depression for tonic neurons
%     m(j+1,n1+1:end) = 0;
%     n(j+1,n1+1:end) = 0;
    
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
        %if Options.doplot1 || Options.doplot2
            raster(j,a) = 1;
        %end
        
    end
    
    spikes(j+1) = spikes(j)+ss/N+dt*(-spikes(j)/tauavg);
    
    if Options.trackstatevariablesmeans
        meany(j+1) = mean(y(j+1,:));
        meansi(j+1) = mean(si(j+1,:));
    end
    
    
    %%% giffing %%%%
    if mod(j-1,20*gifresolution) == 0 && dogifplot
        if j > 20*3000 % over 2 seconds
            [P,P2,P3] = updategif(theta(j,:),pmin,pmax,Isummed,gf,gffilename,j,P,P2,P3,raster,spikes,t,vin,istate);
        end
    end
    
end

%% plotting
disp('Plotting....')

if doplot1
    plot1(Parameters,Options,t,y,m,n,si,spikes,raster,rr,Ihistory,sihistory)
end

close(wb)

%% save last conditions
lc.theta = theta(end,:);
lc.y = y(end,:);
lc.m = m(end,:);
lc.n = n(end,:);
lc.si = si(end,:);
lc.N = N;

save('lastconditions.mat','lc')

%% generate outputs

if Options.trackstatevariablesmeans
    outputs.meany = meany;
    outputs.meansi = meansi;
end

outputs.spikes = spikes;
outputs.spikeshappy = spikes/tauavg*1000;
outputs.A = A;
outputs.killlist = killlist;
outputs.raster = raster;
outputs.t = t;
outputs.thappy = t/1000;
outputs.theta = theta;

disp('Done with simulate_v1')


%% ODE functions

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

    function noiseout = noise(deltat,NumNeurons,sigma,getnoisy)
        
        if any(getnoisy==false)
            noiseout = sigma * sqrt(deltat) * randn(NumNeurons,1).*getnoisy;
        else
            noiseout = sigma * sqrt(deltat) * randn(NumNeurons,1);
        end
        
    end


end