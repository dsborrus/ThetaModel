function DBPlot_v3(dt,tmax,t,y,m,n,sij,spikes,raster,vin,istate)

disp('Entering plotting function')

doraster = 1; % 2 colord, 1 no color, 0 no raster

D = vin(1);
D = vin(1);
isig1 = vin(2);
iu1 = vin(3);
isig2 = vin(4);
iu2 = vin(5);
N = vin(6);
prob = vin(7);
tauavg = vin(8);
mgain = vin(9);
taum = vin(10);
nrise = vin(11);
taun = vin(12);
sigain = vin(13);
tausi = vin(14);
tautheta=vin(15);
noisesigma=vin(16);
conmat = vin(17);
doablate = vin(18);
t2ablat = vin(19);
numN2kill_wonepulse = vin(20);


set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultlinemarkersize',10);

nplots = 4;

fig = figure('Position',[500 500 1500 700]);
 
% % % plot 1 - network spike trace

ax1 = subplot(nplots,1,1);
hold on
plot(t/1000,(spikes./tauavg)*1000);
% plot std and 2std
astd = std((spikes./tauavg)*1000);
tmean = mean((spikes./tauavg)*1000);
%plot(t/1000,zeros(length(t),1)+tmean+2*astd,'color',[0.6 0.6 0.6])
%plot(t/1000,zeros(length(t),1)+tmean,'color','k')
%plot(t/1000,zeros(length(t),1)+2,'color',[0.6 0.6 0.6])
title('Network Activity')

disp('First plot done, starting raster calculations')

% % % plot 2 - raster

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

    mxspks = max(max(laggyraster))

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

title('Raster Plot')

% % % % % % % % % % % % % % % % % % 
% synaptic depression

if 1
    ax3 = subplot(nplots,1,3);
    plot(t/1000,mean(y,2))
    title(['Mean synaptic depression (y) of network'])
else
    ax3 = subplot(nplots,1,3);
    plot(t/1000,y(:,rr))
    title(['One random neuron`s synaptic depression (y) (n=' mat2str(rr) ')'])
end

linkaxes([ax1 ax2 ax3],'x')

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
       '. t2ablat = ' mat2str(t2ablat) '. numN2kill_wonepulse = ' mat2str(numN2kill_wonepulse)];
annotation('textbox',[.1 .03 .1 .1],'String',str,'FitBoxToText','on');

axis off

% lala = length(dir('output/*.png'))+1;

% saveas(fig,['output/dataset' mat2str(lala) '.png'])

% save(['output/dataset' mat2str(lala) '.mat'])
% 
% figure
% plot(t/1000,sij(:,rr))
% title(['One random neuron`s synaptic current output (s) (n=' mat2str(rr) ')'])
