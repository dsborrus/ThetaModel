% Sync Theta Model - Version 5 DB version
% Code written by Greg, commented by Dan, further edits by Dan
% 
% System of linked theta neurons, with two populations of different
% intrinsic frequency
%
% Function is
% d theta/ d time = 1 - cos(theta) + [1 + cos(theta)] * I 
% hint: positive I is instrinicly rhymic

clc
clear all
close all

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

% % % User Params % % %

n1 = 100;   % number of neurons in the first population
n2 = 0;     % number of neurons in the second population
dt = 0.1;   % time step
tnum = 5e4; % the nummber of time steps
iu1 = 0.01;  % mean I parameter for first population
isig1 = 0.001;  % std of I parameter for first population
iu2 = 2.0;  % mean I parameter for second population
isig2=0.0;    % std of I parameter for second population
prob = 0.9; % E-R graph, prob is prob of connection.
D=0.1;      % Strength of networkness
tauavg=1;   % Relaxation of network excitement

% % % Script Settings % % %

DoGif = 0; % Should we show animation? And save the gif?
nplot=5;   % How long to wait before plotting gif
DoDBPlot = 1; % Should we make the graph dan made?
Gradeit = 1; % Use fft to grade the output (NetVolAnalysis)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

n=n1+n2;    % total pop number
pmin = -pi; % domain min
pmax = pi;  % domain max
t=dt*(1:tnum); % sneaky way to make the time array, with correct t steps
phase = pmin:0.01:pmax; 
theta = zeros(tnum,n); % initialize the ODE out put (theta array)
spikes = NaN*ones(tnum,1); spikes(1)=0; % initialize the spike array
rplot = NaN*ones(tnum,1); rplot(1)=0;   % initialize the raster plot?
I = [ iu1+isig1*randn(1,n1) iu2+isig1*randn(1,n2) ];  % initialize I vector
                                                      % recall, it should 
                                                      % be length of n1+n2
delta = D/n; % strength of connections (total excitability of network 
             % divided by the number of connections)
A=(rand(n,n)<prob); % creates connectivity matrix, damn that's savy

% makes sure no neuron is connected to itself. Again, quite savy
for m=1:n, A(m,m)=0; end

% set intialize values for theta, it appears they can start outside -pi to 
% pi. I will fix that.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old Greg stuff:
%theta(1,:) = (pmax+pmin)/2+(pmax-pmin)*randn(1,n);
%theta(1,:) = 0.1*randn(1,n);
% new Dan stuff:
theta(1,:) = pmin + (pmax-pmin).*rand(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set name for gif output
basename = [ 'fig ' datestr(now) ];
% set name
filename = [ basename '.gif' ];
% important again for gif making, need to know if we are on first gif image
firstframe = 1;

% Begin simulation loop
for j = 1:tnum-1
    
    % This is the ODE
    theta(j+1,:) = theta(j,:)+dt*(1-cos(theta(j,:))+(1+cos(theta(j,:))).*I );
    % 
    e=1; ss=0;
    while any(e)
    % check if any neuron is above pi ( should be pmax )
    e = (theta(j+1,:)>pi); a=find(e);
    % reset that neuron back 2pi
    theta(j+1,a)=theta(j+1,a)-2*pi;
    % calculate the number of pulses the neuron receives
    s=e*A;
    % Add the pulse
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % EDITS, Greg was adding delta*s*the phse response curve
    % I will chage it to only add delta*s
    % Greg's version:
    %theta(j+1,:) = theta(j+1,:)+delta*s.*(1+cos(theta(j+1,:)));
    % Dan's version (new):
    theta(j+1,:) = theta(j+1,:)+delta*s;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    ss = ss+sum(s);
    end

    spikes(j+1) = spikes(j)+dt*(ss/n-spikes(j)/tauavg);
    rplot(j+1) = 1/n*sum(exp(i*theta(j+1,:)));
    
    if j-nplot>1 && DoGif
        PlotAndGifSave
    end
end

if DoGif
    print([ basename '.pdf'],'-dpdf')
end

if DoDBPlot
    DBPlot
end

if Gradeit
    
    % wait 100 seconds
    eqspikes = spikes(100/dt+1:end);
    NetVolAnalysis(eqspikes,dt);
end

% fid = fopen(['figdata.dat'],'w');
% for k=1:length(t)
%     shift=100;
%     if t(k)>=shift
%     fprintf(fid,'%6.4f  ',[t(k) theta(k,:) w(k,:)]);
%     fprintf(fid,'\n');
%     end
% end
% fclose(fid);


