function plot1(Parameters,Options,t,y,m,n,sij,spikes,raster,rr,Ihistory,sijhistory)

D = Parameters.D;
isig1 = Parameters.isig1;
iu1 = Parameters.iu1;
isig2 = Parameters.isig2;
iu2 = Parameters.iu2;
N = Parameters.n1 + Parameters.n2;
prob = Parameters.prob;
tauavg = Parameters.tauavg;
mgain = Parameters.mgain;
taum = Parameters.taum;
nrise = Parameters.nrise;
taun = Parameters.taun;
sigain = Parameters.sigain;
tausi = Parameters.tausi;
tautheta=Parameters.tautheta;
noisesigma=Parameters.noisesig;
dt = Parameters.dt;
tmax = Parameters.tmax;
istate = Parameters.istate;
conmat = Options.conmat;
doablate = Options.Doablate;
t2ablat = Parameters.t2ablat;
N2k_w1p = Parameters.N2k_w1p;
doraster = Options.doraster;

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

nplots = 9;

tstart = dt;
tend = tmax;

tw = tstart/dt:tend/dt;

fig = figure('Position',[800 500 1500 1200]);

ax1 = subplot(nplots,1,1);
hold on
plot(t(tw)/1000,(spikes(tw)./tauavg)*1000)
title('Network Activity')

%%%%%%%%%%%%%%%%% plot 2 - raster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = subplot(nplots,1,2); hold on;


if doraster == 2
    win = 700 /dt; % window to use to capture number of spikes in that window (/dt should be there for conversions to framerate)
    laggyraster = zeros(length(t),N);

    for i = 2:length(t)
        if i<=win
            laggyraster(i,:) = sum(raster(1:i,:),'omitnan');
        else
            laggyraster(i,:) = sum(raster(i-win:i,:),'omitnan');
        end
    end

    laggyraster = laggyraster .* raster;

    mxspks = max(max(laggyraster));

    disp('Done calculating raster, now plotting it')
    % % % % % % % %      % % % % % % %
    % first colormap is just analog color scale
    %colormap = [linspace(0,1,mxspks)' zeros(mxspks,1) zeros(mxspks,1)];
    % second colormap is discrete kinda random colors, will break if too many
    % spikes
    colormap = [0 0 0;      % black         1
                0 0 0.5     % dark blue     2
                0 0.5 0     % dark green    3
                0 0 1       % blue          4
                0 1 0       % green         5
                0 1 1       % cyan          6 
                1 0 1       % magenta       7
                1 0 0       % red           8
                1 0.5 0     % orange        9
                1 1 0       % yellow        10
                1 0.8 0.8   % pink          11
                ];


    for c = 1:mxspks

        templaggyraster = laggyraster==c;

        for i = 1:N
            plot(t/1000,templaggyraster(:,i)*i,'.','color',colormap(c,:)); hold on;
        end

    end

elseif doraster == 1

    for i = 1:N
        plot(t/1000,raster(:,i)*i,'.k'); hold on;
    end
    
end
    
ylim([1 N])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title('Raster Plot')

if 1
ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,mean(y(tw,:),2))
title(['Mean synaptic depression (y) (n=' mat2str(rr) ')'])
else
ax3 = subplot(nplots,1,3);
plot(t(tw)/1000,y(tw,rr))
title(['One random neuron`s synaptic depression (y) (n=' mat2str(rr) ')'])
end

ax4 = subplot(nplots,1,4);
plot(t(tw)/1000,sij(tw,rr))
title(['One random neuron`s synaptic current output (s) (n=' mat2str(rr) ')'])

ax5 = subplot(nplots,1,5);
plot(t(tw)/1000,m(tw,rr))
title(['One random neuron`s (m) (n=' mat2str(rr) ')'])

ax6 = subplot(nplots,1,6);
plot(t(tw)/1000,n(tw,rr))
title(['One random neuron`s (n) (n=' mat2str(rr) ')'])

ax7 = subplot(nplots,1,7);
plot(t(tw)/1000,sijhistory(tw))
title(['One random neuron`s synaptic current input (n=' mat2str(rr) ')'])

ax8 = subplot(nplots,1,8);
plot(t(tw)/1000,Ihistory(tw))
title(['One random neuron`s total applied current (I+delta*syncurrent*syndepression) (n=' mat2str(rr) ')'])

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'x')

subplot(nplots,1,nplots);
hold on

str = ['D = ' mat2str(D) '. isig1 = ' mat2str(isig1) '. iu1 = ' mat2str(iu1) ...
    '. isig2 = ' mat2str(isig2) '. iu2 = ' mat2str(iu2) ...
    '. N = ' mat2str(N) '. prob = ' mat2str(prob) '. tauavg = ' mat2str(tauavg) ...
    '. mgain = ' mat2str(mgain) '. taum = ' mat2str(taum) ...
    '. sigain = ' mat2str(sigain) '. tausi = ' mat2str(tausi)];
annotation('textbox',[.1 .08 .1 .1],'String',str,'FitBoxToText','on');

str = ['nrise = ' mat2str(nrise) '. taun = ' mat2str(taun) ...
       '. tautheta = ' mat2str(tautheta) '. noisesigma = ' mat2str(noisesigma)...
       '. istate = ' mat2str(istate) '. conmat = ' mat2str(conmat) '. doablate = ' mat2str(doablate)...
       '. t2ablat = ' mat2str(t2ablat) '. numN2kill_wonepulse = ' mat2str(N2k_w1p)];
annotation('textbox',[.1 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

%saveas(fig,['dataset' mat2str(length(dir('*.png'))+1) '.png'])
