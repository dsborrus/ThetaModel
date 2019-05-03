% Sync Theta Model with synaptic depression - Version 3 DB version
% Orignal ode written by Greg, commented by Dan, further edits by Dan
% Original ode was simple theta model This now has synaptic depresseion
%
% % % % Update History % % % %
function [yss,sps] = synctheta_function_noy(D,ydrop,tmax,tauy,istate,bumpit,n1,n2,iu1,iu2,isig1,isig2,plotit)
% This version has been remade into a function
% This function should return the steady state value of y, if y is not
% connected to any part of the function 
% i.e. the pulse is just deltatheta. NOT deltatheta * y
% the value ydrop should be an important parameter in making the desired
% shape of the curve
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

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

% % % User Params % % %
% 
% n1 = 180;   % number of neurons in the first population
% n2 = 20;     % number of neurons in the second population
dt = 0.01;   % time step
%tmax = 2e3;     % maximum time of simulation
% iu1 = -0.01;  % mean I parameter for first population
% isig1 = 0.00;  % std of I parameter for first population
% iu2 = 0.01;  % mean I parameter for second population# 
% isig2 = 0.000;    % std of I parameter for second population
prob = 0.5; % E-R graph, prob is prob of connection.
%D = 3.1;      % Strength of networkness
tauavg=1;   % Relaxation of network excitement

%ydrop = .2; % How much of an effect firing has on synaptic depression
            % (should be between 0 and 1)!!!
%tauy  =  6; % Char time for return to ss for y (synap depress)     

dur = 30;
           
tstart_an = tmax/2;

tend_an = tmax-1;

% % % Script Settings % % %

DoDBPlot = plotit; % Should we make the graph dan made?

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

raster = NaN(tnum,n); % initialize the raster plot

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
        theta(1,:) = pmin+randn(1,n)*.01;
        y(1,:) = (randn(1,n)*.01)+0.1;
    case 2 %high state
        theta(1,:) = pmax.*rand(1,n);
        y(1,:) = (randn(1,n)*.01)+.97;
end

% initialize spike counter
% spikes per frame
spf = zeros(tnum,1);

% Begin simulation loop
for j = 1:tnum-1
    
    if bumpit
        if isittime(t(j),dur)
            Ibump = rand(1,n)*0.5;
            % update theta, we are bumping
            theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),I+Ibump);
        else
            % update theta, we aren't bumping anymore 
            theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),I);
        end
    else
        % update theta, we will never bump
        theta(j+1,:) = theta(j,:) + dt * thetaODE(theta(j,:),I);
    end
    
    % Calculate ODEs next step (Euler's method)
    y(j+1,:)     = y(j,:)     + dt * yODE(y(j,:),tauy);
    
    
    e=1; ss=0;
    while any(e)
    % check if any neuron is above pi ( should be pmax )
    e = (theta(j+1,:)>pi); a=find(e);
    % reset that neuron back 2pi
    theta(j+1,a)=theta(j+1,a)-2*pi;
    % augment the synpatic depression term for next y
    y(j+1,a) = y(j+1,a) - ydrop;
    
    % I grab rows from A and multiply them by 
    % associated synaptic depression. Then sum them along columns
    % We only grab the rows and columns we need (to speed up code)
    %s = sum(A(e,:).*y(j,e)',1);
    
    % uncouple pulse from y
    s = sum(A(e,:),1);
    
    
    % Calculate and add the pulse
    theta(j+1,:) = theta(j+1,:)+delta*s;
    
    % Sum up number of spikes in this step
    ss = ss+length(a);
    
    % catalouge for raster plot
    if DoDBPlot
        raster(j,a) = 1;
    end
    
    % little check for sanity %
    if any(delta*s>2*pi)
        warning('I think this means we pushed a neuron up 2pi, it missed a phase')
    end
    
    end

    spikes(j+1) = spikes(j)+dt*(ss/n-spikes(j)/tauavg);
    %rplot(j+1) = 1/n*sum(exp(i*theta(j+1,:)));
    spf(j) = ss;
    
end

toc

%% DB plot %%

if DoDBPlot
    
    tstart_pl = dt;
    tend_pl = tmax;


if tstart_pl>=tmax
    tstart_pl = dt;
end

if tend_pl>tmax
    tend_pl = tmax;
end

tw_pl = tstart_pl/dt:tend_pl/dt;

figure('Position',[800 500 1300 800])

ax1 = subplot(4,1,1);
hold on
plot(t(tw_pl),spikes(tw_pl))
title('Network Activity')


ax2 = subplot(4,1,2);
temp = dt:dt:tmax;
for i = 1:n
    plot(temp(tw_pl),raster((tw_pl),i)*i,'k.'); hold on;
end
title('Raster Plot')

ax3 = subplot(4,1,3);
rr = randi(n);
plot(t(tw_pl),y(tw_pl,rr))
title(['One random neuron`s synaptic depression (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3],'x')

subplot(4,1,4)
hold on
str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. n = ' mat2str(n) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. tauy = ' mat2str(tauy) '. ydrop = ' mat2str(ydrop)];
annotation('textbox',[.1 .1 .1 .1],'String',str,'FitBoxToText','on');

axis off

end


%% Analysis %%

tw_an = tstart_an/dt:tend_an/dt;

spf_insidewindow = spf(tw_an);

sps = sum(spf_insidewindow)/(length(spf_insidewindow)*dt);

% % %

% how to calculate y_ss... take average value for last 10 seconds?
% of all neurons
y_an_t = 10;

yss = mean(mean(y(end-(y_an_t/dt):end,:)));

%% Other functions

function dtheta = thetaODE(theta,I)
    dtheta = 1-cos(theta) + (1+cos(theta)).*I;
end

function dy = yODE(y,tauy)
    dy = (1 - y)/tauy;
end

function YorN = isittime(t,dur)
    
    time1 = 100;
    time2 = 200;
    time3 = 300;
    
    if          t>time1 && t<time1+dur
        YorN = 1;
    elseif      t>time2 && t<time2+dur
        YorN = 1;
    elseif      t>time3 && t<time3+dur
        YorN = 1;
    else
        YorN = 0;
    end
end

end
