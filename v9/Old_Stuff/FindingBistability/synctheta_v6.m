function synctheta_v6(params,h)
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

n1 = params.n1;   % number of neurons in the first population
n2 = params.n2;     % number of neurons in the second population
dt = params.dt;   % time step
tmax = params.tmax; % maximum time of simulation
iu1 = params.iu1;  % mean I parameter for first population
isig1 = params.isig1;  % std of I parameter for first population
iu2 = params.iu2;  % mean I parameter for second population#
isig2 = params.isig2;    % std of I parameter for second population
prob = params.prob; % E-R graph, prob is prob of connection.
D = params.MaxD;      % Strength of networkness
tauavg=params.tauavg;   % Relaxation of network excitement

sigain = params.sigain;
tausi = params.tausi;

ydrop = params.ydrop; % How much of an affect firing has on synaptic depression
% % (should be between 0 and 1)!!!
tauy = params.tauy; % Char time for return to ss for y (synap depress)

silence = params.silence;
bumpit = params.bumpit;
istate = params.istate;

% plot stuff

subplot(5,3,[7 8]); hold on;
pl1 = plot(0,0);
pl2 = plot(0,0);
pl3 = plot(0,0);
pl4 = plot(0,0);
% pl5 = plot(0,0);
ylim([0 1]);
axis off
subplot(5,3,[10 11]); hold on;
ylim([1 n1+n2]);
wind = 100;

temp = length(dir('*.gif'));
filename = (['AnimatedGIF' mat2str(temp) '.gif']);

% % % Script Settings % % %

DoDBPlot = 1; % Should we make the graph dan made?
dogif = 0;

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
Ihistory = NaN(tnum,2);
sihistory = NaN(tnum,2);
%ysstarray = zeros(tnum,1);
%spstarray = zeros(tnum,1);
% random neuron to track
rr1 = randi(n1);
rr2 = randi(n1);
yssprev = [NaN NaN NaN];
spsprev = [NaN NaN NaN];
%yssn2_prev = [NaN NaN NaN];
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
    Isummed = I + (delta*y(j,:).*si(j,:)) * A; 
    %Isummed = I + (delta*sij(j,:)) * A; 
    
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
    Ihistory(j,:) = Isummed([rr1 rr2]);
    sihistory(j,1) = si(j,:) * A(:,rr1);
    sihistory(j,2) = si(j,:) * A(:,rr2);

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
    
    spikes(j+1) = spikes(j)+dt*(ss/n-spikes(j)/tauavg);
    
    spf(j) = ss;
    

    
if t(j) >= 2e3
if mod(t(j),1000) == 0
    
    %rc = mod(rc,3)+1;
    
    yss = mean(mean(y(j-1000/dt:j,:)));
    yssprev(end+1) = yss;
    yssprev(1) = [];
    
    yssall = mean(y(j-1000/dt:j,:));
    
%     yssn2 = mean(mean(y(j-1000/dt:j,n1+1:end)));
%     yssn2_prev(end+1) = yssn2;
%     yssn2_prev(1) = [];
    
    sps = sum(spf(j-1000/dt:j));
    spsprev(end+1) = sps;
    spsprev(1) = [];
   
    
    subplot(5,3,[1 6]); hold on;
    delete(pl1);
    delete(pl3);
    pl1 = plot(yss,sps,'-','color',[0 0 1],'HandleVisibility','off');
    pl3 = plot(yssprev,spsprev,'o','color',[0 0 1],'HandleVisibility','off');
    
    delete(pl4)
    pl4 = plot(yssall,sps','.','color','k','HandleVisibility','off');
    
%     delete(pl4);
%     delete(pl5);
%     pl4 = plot(yssn2,sps,'-','color',[0 0.5 0]);
%     pl5 = plot(yssn2_prev,spsprev,'o','color',[0 0.5 0]);
    %legend('Low init state','High init state + mid simulation bump','y_{ss}','','');

        
     subplot(5,3,[7 8]); hold on;
     delete(pl2)
     pl2 = plot(j/dt-wind:j/dt,spikes((j-wind):j),'b');
    
    subplot(5,3,[10 11]); hold off;
    for k = 1:(n1+n2)
        plot(j/dt-wind:j/dt,raster(j-wind:j,k)*k,'k.'); hold on;
    end
    
    subplot(5,3,[13 15])
    
    delete(findall(gcf,'Tag','somethingUnique'))
    str = ['t = ' mat2str(t(j))];
    annotation('textbox',[.7 .4 .1 .1],'String',str,'FitBoxToText','on'...
        ,'Tag' , 'somethingUnique');
    
    drawnow
    

    %pause(1)
    
        if dogif
        % capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write the GIF File
        if t(j) == 2e3
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append')
        end
        end
end  
end
    
    
    
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
    
    DBPlot_v2(dt,tmax,t,n,y,si,spikes,raster,rr1,Ihistory(:,1),sihistory(:,2),vin)
    DBPlot_v2(dt,tmax,t,n,y,si,spikes,raster,rr2,Ihistory(:,1),sihistory(:,2),vin)
end

%% Analysis %%

% time for analysis window
tstart_an = tmax/2;
tend_an = tmax-1;

tw_an = tstart_an/dt:tend_an/dt;

spf_insidewindow = spf(tw_an);

sps = sum(spf_insidewindow)/(length(spf_insidewindow)*dt/1000);


%%%%

yss = mean(mean(y(:,end-end/2)));

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