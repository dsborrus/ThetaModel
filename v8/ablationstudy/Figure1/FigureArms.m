classdef FigureArms < handle
    %FIGUREPARTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Parameters
        Options
        str
        A
        Aexists
        
        simplet
        ablatet
        outsimple
        outablate
        
        
    end
    
    methods
        function obj = FigureArms(parameters,options,str)
            %FIGUREPARTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Parameters = parameters;
            obj.Options = options;
            obj.str = str;
            obj.Aexists = 0;
        end
        
        function MakeGraph(obj)
            
            addpath('../..')
            
            obj.A = MakeNetwork(obj.Parameters,obj.Options);
            
            obj.Aexists = 1;
            
        end
        
        function AnalyseGraph(obj,extraplot_in,extraplot_out)
            %MakeConnectGraph Make Connectivity Matrix
            %   This method creates and saves the connectivity matrix of a
            %   given conmat input
            %   conmat = 1 - er
            %   conmat = 2 - sw
            %   conmat = 3 - sf
            %   conmat = 5 - ke
            
            
            % e4xtraplot 0 = no extra lines
            %extra plot 1 = possoin dist
            % extra plot 2 = power
            
            % Connectivity Matrix Figure
            
            if obj.Aexists==0
                error('must build network first (run MakeGraph)')
            end
            
            N = obj.Parameters.n1+obj.Parameters.n2;
            if ~exist('extraplot_in','var')
                extraplot_in = 0;
            end
            
            if ~exist('extraplot_out','var')
                extraplot_out = 0;
            end
            
            
            f = figure; hold on
            spy(obj.A);

            ax = gca;

            ax.FontSize = 16;
            ax.XTick = [1 100];
            ax.YTick = [1 100];

            xlabel('')

            saveas(f,['AdjacencyMatrixFigures/' obj.str '.png'])
            
            close all
            
            % 2 - make kin figures
            
            kin = sum(obj.A);
            kout = sum(obj.A,2);
            
            % convert into probability degree distribution
            pkin = histcounts(kin,max(kin)-min(kin)+1)/(N);
            pkout = histcounts(kout,max(kout)-min(kout)+1)/(N);
            
            % also find the mean in/out degree
            meankin = mean(kin);
            meankout = mean(kout);
            
%             [in_p,in_k_pre] = histcounts(kin);
%             [out_p,out_k_pre] = histcounts(kout);
%             
%             in_k = movmean(in_k_pre,2);
%             in_k = in_k(2:end);
%             out_k = movmean(out_k_pre,2);
%             out_k = out_k(2:end);
%             
%             indata = [in_k' in_p'];
%             outdata = [out_k' out_p'];

            indata = [(min(kin):max(kin))' pkin'];
            outdata = [(min(kout):max(kout))' pkout'];
            
            poif = @(k,lambda) ( (lambda.^k)*(exp(-lambda)) )./factorial(k);
            powf = @(k,kmin,gamma) (gamma-1)*kmin^(gamma-1)*k.^(-gamma);
            
            
            f = figure;
            if extraplot_in == 0
                plot((min(kin):max(kin)), pkin, 'o','markersize',15,'markerfacecolor','b'); hold on;
            elseif extraplot_in == 1
                plot((min(kin):max(kin)), pkin, 'o','markersize',15,'markerfacecolor','b'); hold on;
                k = (0:max(kin)+5);
                plot(k,poif(k,N*obj.Parameters.er_prob),'linewidth',3)
            elseif extraplot_in == 2
                loglog((min(kin):max(kin)), pkin, 'o','markersize',15,'markerfacecolor','b'); hold on;
                k = (min(kin)-5:max(kin)+5);
                k = min(xlim):max(xlim);
                loglog(k,powf(k,min(kin),3),'linewidth',3)
            end
            theCA = get(gca);
            ylim([theCA.YLim(1) theCA.YLim(2)])
            if extraplot_in == 2
                plot([meankin meankin],[theCA.YLim(1) 10^(log10(theCA.YLim(1)) + (0.2*( log10(theCA.YLim(2)) - log10(theCA.YLim(1)) ) ))],'r','linewidth',2)
            else
                plot([meankin meankin],[theCA.YLim(1) (theCA.YLim(2)-theCA.YLim(1))*.20],'r','linewidth',2)
            end
            set(gca,'fontsize',23)
            saveas(f,['KinData/' obj.str '.png'])
            close all
            
                        
            f = figure;
            if extraplot_out == 0
                plot((min(kout):max(kout)), pkout, 'o','markersize',15,'markerfacecolor','b'); hold on;
            elseif extraplot_out == 1
                plot((min(kout):max(kout)), pkout, 'o','markersize',15,'markerfacecolor','b'); hold on;
                k = (0:max(kout)+5);
                plot(k,poif(k,N*obj.Parameters.er_prob),'linewidth',3)
            elseif extraplot_out == 2
                loglog((min(kout):max(kout)), pkout, 'o','markersize',15,'markerfacecolor','b'); hold on;
                k = (min(kout)-5:max(kout)+5);
                k = min(xlim):max(xlim);
                loglog(k,powf(k,min(kout),3),'linewidth',3)
            end
            theCA = get(gca);
            ylim([theCA.YLim(1) theCA.YLim(2)])
            if extraplot_out == 2
                plot([meankout meankout],[theCA.YLim(1) 10^(log10(theCA.YLim(1)) + (0.2*( log10(theCA.YLim(2)) - log10(theCA.YLim(1)) ) ))],'r','linewidth',2)
            else
                plot([meankout meankout],[theCA.YLim(1) (theCA.YLim(2)-theCA.YLim(1))*.20],'r','linewidth',2)
            end
            set(gca,'fontsize',23)
            saveas(f,['KoutData/' obj.str '.png'])
            close all

%             save(['KinData/' obj.str '.dat'],'indata','-ascii')
%             save(['KoutData/' obj.str '.dat'],'outdata','-ascii')
            
        end
        
        function simplsimNetwork(obj,loadoldsim)
            % network is built, lets just run with NO ablation
           
            % loadoldsim determines if we should load an old sim, or run a
            % new one. We can save time by loading an old one. But if we
            % changes params, we dont want to use and old run :/
            % loadoldsim == 1 means we do load an old run
            % loadoldsim == 0 (default) means we rerun sim
           
            % overwriting tmax!!!!!! Just because :/
            obj.Parameters.tmax = 60000;
            obj.Options.Doablate    = 0;
            obj.Options.trackstatevariablesmeans=0;

            if obj.Aexists==0
                error('must build network first (run MakeGraph)')
            end

            % assigning matrix as a parameter, so that we dont generate a
            % new network
            obj.Parameters.A = obj.A;

            addpath('../..')

            % simulate
            if ~exist('loadoldsim')
                loadoldsim = 0;
            end
            % options to load old sims to save time
            if loadoldsim == 0
                out = simulate_v1(obj.Parameters,obj.Options);
            else % load old sim
                if isempty(obj.outsimple)
                    error('Dude, you gotta run the sim once to have one saved')
                else
                    out = obj.outsimple;
                end
            end
            
            obj.simplet = obj.Parameters.dt:obj.Parameters.dt:obj.Parameters.tmax;
            t = obj.simplet;
                

            % plot
            tstart = 1500;
            tend = obj.Parameters.tmax;
            dt = obj.Parameters.dt;
            tauavg = obj.Parameters.tauavg;

            tw = tstart/dt:tend/dt;

            f = figure; hold on;
            plot(t(tw)/1000,(out.spikes(tw)./tauavg)*1000,'linewidth',2)
            xlim([t(tw(1))/1000 t(tw(end))/1000])
            ylim([0 20])
            %title('Network Activity (intact)')
            xlabel('Time (s)')
            ylabel('Spikes/neuron/second')
            set(gca,'fontsize',23)
            saveas(f,['IntactTraces/' obj.str '.png'])
            close all

            obj.outsimple = out;            
            
            
        end    
        
        function simublateNetwork(obj,loadoldsim)
            % the network is built, time to run it and ablate it

            
            % ablation params %
                       % loadoldsim determines if we should load an old sim, or run a
            % new one. We can save time by loading an old one. But if we
            % changes params, we dont want to use and old run :/
            % loadoldsim == 1 means we do load an old run
            % loadoldsim == 0 (default) means we rerun sim
           
            % overwriting tmax!!!!!! Just because :/
            obj.Parameters.tmax = 1500000; %(in ms)
            %obj.Parameters.tmax = 20000; %(in ms)
            
            % overwriting ablation params as well!!!
            obj.Parameters.t2ablat  = 300;     % When to ablate neurons (every XXX seconds)
            obj.Parameters.N2k_w1p  = 10;       % number of neurons to kill with one pulse
            
            % remember to ablate
            obj.Options.Doablate = 1;
            
            % track state variables
            obj.Options.trackstatevariablesmeans=1;

            if obj.Aexists==0
                error('must build network first (run MakeGraph)')
            end

            % assigning matrix as a parameter, so that we dont generate a
            % new network
            obj.Parameters.A = obj.A;

            addpath('../..')
            
            % simulate
            if ~exist('loadoldsim')
                loadoldsim = 0;
            end
            % options to load old sims to save time
            if loadoldsim == 0
                out = simulate_v1(obj.Parameters,obj.Options);
            else % load old sim
                if isempty(obj.outablate)
                    error('Dude, you gotta run the sim once to have one saved')
                else
                    out = obj.outablate;
                end
            end
            
            obj.ablatet = obj.Parameters.dt:obj.Parameters.dt:obj.Parameters.tmax;
            t = obj.ablatet;
                

            % plot
            tstart = 1000;
            tend = obj.Parameters.tmax;
            dt = obj.Parameters.dt;
            tauavg = obj.Parameters.tauavg;

            tw = tstart/dt:tend/dt;

            f = figure; hold on;
            plot(t(tw)/1000,(out.spikes(tw)./tauavg)*1000)
            xlim([t(tw(1))/1000 t(tw(end))/1000])
            ylim([0 20])
            plot(obj.Parameters.t2ablat:obj.Parameters.t2ablat:obj.Parameters.tmax/1000,18,'r.','markersize',30)
            %title('Network Activity (intact)')
            xlabel('Time (s)')
            ylabel('Spikes/neuron/second')
            set(gca,'fontsize',23)
            saveas(f,['AblateTraces/' obj.str '.png'])
            close all

            obj.outablate = out;    
           
        end
        
        function plotchunksofablation(obj,chunknumber,window,windowfrontbuffer)
            % this function plots zoomed in chunks of ablation traces
            % chunknumber is the chunk that you want to zoom in on 
            % chunknumbers start at 0 (no ablations)
            % window is the size of the zoom in (length of zoom in seconds)
            % windowfrontbuffer is the amount of room to leave infront of
            % the window (the buffer, ya dunce!!)
            
            if isempty(obj.outablate)
                error('You need to simulate the ablation experiment first')
            end
            if window > obj.Parameters.t2ablat
                error('Your window (in s) is longer than the total ablation window')
            end
            if window+windowfrontbuffer > obj.Parameters.t2ablat
                error('Your window and/or windowfrontbuffer are too large (larger than t2ablat')
            end
            
            t2ablat = obj.Parameters.t2ablat;     % When to ablate neurons (every XXX seconds)
            dt = obj.Parameters.dt;
            tauavg = obj.Parameters.tauavg;
            
            tstart_s = (chunknumber*t2ablat) + windowfrontbuffer;
            tstart_j = tstart_s*1000/dt;
            
            tend_s = (chunknumber*t2ablat) + windowfrontbuffer + window;
            tend_j = tend_s*1000/dt;
            
            tw = tstart_j:tend_j;
            
            f = figure; hold on;
            plot(obj.ablatet(tw)/1000,(obj.outablate.spikes(tw)./tauavg)*1000,'linewidth',2)
            xlim([obj.ablatet(tw(1))/1000 obj.ablatet(tw(end))/1000])
            %ylim([0 20])
            %plot(obj.Parameters.t2ablat:obj.Parameters.t2ablat:obj.Parameters.tmax/1000,18,'r.','markersize',30)
            %title('Network Activity (intact)')
            xlabel('Time (s)')
            ylabel('Spikes/neuron/second')
            set(gca,'fontsize',23)
            
            % hardcoding in what folder to save to
            switch chunknumber
                case 1
                    ylim([0 15])
                    saveas(f,['10percent_ablated/' obj.str '.png'])
                case 2
                    ylim([0 15])
                    saveas(f,['20percent_ablated/' obj.str '.png'])
                case 3
                    ylim([0 15])
                    saveas(f,['30percent_ablated/' obj.str '.png'])
                case 4
                    ylim([0 15])
                    saveas(f,['40percent_ablated/' obj.str '.png'])
            end
                    
            close all
            
            
        end
        
        function AblationRhythmSummary(obj,Fs)
            % This function creates a figure comparing all the summary
            % statistics of rhythm breakdown
            
            % obj is of course the first F1
            % but then Fs is an column array of all the networks to compare
            % it to (including F1) ([F1; F2; F3; ... F4])
            
            if isempty(obj.outablate)
                error('You need to simulate the ablation experiment first')
            end
            
            N2k_w1p = obj.Parameters.N2k_w1p;
            t2ablat = obj.Parameters.t2ablat;
            tmax    = obj.Parameters.tmax;
            dt      = obj.Parameters.dt;
            
            % calculate maximum percentage killed
            MaxAblated = N2k_w1p * (floor(tmax/1000/t2ablat) - 1);
            
            if MaxAblated <= 0 
                disp('Something is wrong with tmax or other ablation params')
                disp('Maybe you ran a normal sim after running an ablation sim?')
                disp('Try replotting ablation sim, without resimming')
                error('See error message above')
            end
                
            Chunks = 0:N2k_w1p:MaxAblated;
            
            NormChunk = Chunks./(obj.Parameters.n1+obj.Parameters.n2);
            
            frontwin = 20; % arbitrarily picking here, a front window of 20 
            % seconds for analysis
           
            f1 = figure; hold on;
            f2 = figure; hold on;
            f3 = figure; hold on;
            
            % loop through simulation
            for i = 1:length(Fs)
                
                F = Fs(i);
                fullspks = (F.outablate.spikes./obj.Parameters.tauavg)*1000;
                
                brstthresh1 = mean(fullspks(frontwin*1000/dt:t2ablat*1000/dt));
                brstthresh2 = mean(fullspks(frontwin*1000/dt:t2ablat*1000/dt))*0.8;
                %brstthresh3 = mean(fullspks(frontwin*1000/dt:t2ablat*1000/dt))+2*std(fullspks(frontwin*1000/dt:t2ablat*1000/dt));
                
                
                meanamps = [];
                meanarea = [];
                meanperiod = [];
                
                varamps = [];
                vararea = [];
                varperiod =[];
                
                % loop through chunks
                for j = 0:length(Chunks)-1
                    
                    tstart_s = (j*t2ablat) + frontwin;
                    tstart_j = tstart_s*1000/dt;

                    tend_s = (j*t2ablat) + t2ablat;
                    tend_j = tend_s*1000/dt;

                    tw = tstart_j:tend_j;
                    
                    spks = fullspks(tw);
                    microt = F.ablatet(tw);
                    
                    amps  = [];
                    area  = [];   
                    burstt = [];
                    
                    up = false;
                    % loop through spike points
                    for k = 2:length(spks)
                        
                        % only triggered when moving above 1st threshold
                        % and down
                        if spks(k) > brstthresh1 && up == false
                            up = true;
                            amps(end+1) = spks(k);
                            burstt(end+1) = k;
                        end
                        
                        % only triggered when above and positive slope
                        if up == true && spks(k) > amps(end)
                            amps(end) = spks(k);
                            burstt(end) = k;
                        end
                        
                        % only triggered when already 
                        if spks(k) < brstthresh2 && up == true
                            up = false;
                            endt = k;
                            area(end+1) = sum(spks(burstt(end):endt))*(dt/1000);
                        end
                            
                            
                        
                        
                        
                    end
                    
                    % testing room %
%                     ftest = figure;
%                     plot(microt/1000,spks); hold on;
%                     plot(microt/1000,brstthresh1*ones(length(microt),1));
%                     plot(microt/1000,brstthresh2*ones(length(microt),1));
%                     plot(20+startt*dt/1000,brstthresh1,'x')
                    %plot(microt/1000,brstthresh3*ones(length(microt),1));
                    
                    meanamps(end+1) = mean(amps);varamps(end+1) = std(amps);
                    meanarea(end+1) = mean(area);vararea(end+1) = std(area);
                    meanperiod(end+1)=mean(diff(burstt*dt/1000));varperiod(end+1)=std(diff(burstt*dt/1000));

                end
                
                figure(f1);
                errorbar(Chunks,meanamps,varamps,'o-','markersize',15);
                figure(f2);
                errorbar(Chunks,meanarea,vararea,'o-','markersize',15);
                figure(f3);
                errorbar(Chunks,meanperiod,varperiod,'o-','markersize',15);
                 
            end
            
            
            % last figure stuff!!!
            figure(f1);
            title('Mean amplitude of burst')
            xlabel('Network Destroyed (in %)')
            ylabel('Amplitude')
            legend('ER','SW','SW noisy','SF van','SF curr','S&S van')
            set(gca,'fontsize',20)
            %set(gcf, 'Position',  [100, 100, 2000, 900])
            saveas(f1,['AblationSummary/amp.png'])
            
            figure(f2);
            title('Mean area of burst')
            xlabel('Network Destroyed (in %)')
            ylabel('Area')
            legend('ER','SW','SW noisy','SF van','SF curr','S&S van')
            set(gca,'fontsize',20)
            %set(gcf, 'Position',  [100, 100, 2000, 900])
            saveas(f2,['AblationSummary/area.png'])
            
            figure(f3);
            title('Mean period of rhythm')
            xlabel('Network Destroyed (in %)')
            ylabel('Period (in s)')
            legend('ER','SW','SW noisy','SF van','SF curr','S&S van')
            set(gca,'fontsize',20)
            %set(gcf, 'Position',  [100, 100, 2000, 900])
            saveas(f3,['AblationSummary/per.png'])
            
            close all
                
            
        end
        
        function DrawStateTraces(obj,chunks,vars2comp,figindex)
            % this function creates the figures showing how the state
            % variables are changing
            
            % chunks is an array (ex: [0 2]) of which chunks to analyse
            % array length must be maxchunks>0
            % vars2comp is a length 2 vector with the variables we wish to
            % compare, input as strings. Possible inputs are f,fdot,s,sdot,y,ydot.
            % figindex lets me and latex know how many comparison there
            % are. For just one row of figures, figindex=1, if we want
            % more, fig index can get adjusted
            
            if obj.Options.Doablate ~= 1
                error(['Make sure you run ablation sim before this one, it' ...
                'can even just be replotted, not re-simulated'])
            end
                
            
            if length(vars2comp)~=2
                error
            end
            
            dt = obj.Parameters.dt;
            tmax = obj.Parameters.tmax;
            t2ablat = obj.Parameters.t2ablat;
            
            vars = zeros(tmax/dt,2);
            
            for i = 1:2
                if strcmp(vars2comp(i),'s')
                    vars(:,i) = obj.outablate.meansi;
                elseif strcmp(vars2comp(i),'sdot')
                    vars(2:end,i) = diff(obj.outablate.meansi)./diff(obj.ablatet*dt)';
                elseif strcmp(vars2comp(i),'y')
                    vars(:,i) = obj.outablate.meany;
                elseif strcmp(vars2comp(i),'ydot')
                    vars(2:end,i) = diff(obj.outablate.meany)./diff(obj.ablatet*dt)';  
                elseif strcmp(vars2comp(i),'f')
                    vars(:,i) = (obj.outablate.spikes./obj.Parameters.tauavg)*1000;
                elseif strcmp(vars2comp(i),'fdot')
                    vars(2:end,i) = diff((obj.outablate.spikes./obj.Parameters.tauavg)*1000)./diff(obj.ablatet*dt)';
                else
                    error('You f-ed up')
                end
            end
            
            frontwindow = t2ablat-45; % in seconds
            
            f = figure; hold on;
%             ylim([0.2 1])
%             xlim([0 0.4])
            xlabel(vars2comp(1))
            ylabel(vars2comp(2))
            set(gca,'fontsize',23');
            
            % just for video v
            colors = [0,0,1;
                      0,0,0;
                      0,0,0;
                      1,0,1];
            % just for video ^
            
            for c = chunks
                t1pre = (c*t2ablat) + frontwindow;
                t1 = t1pre*1000/dt;
                t2pre = ((c+1)*t2ablat);
                t2 = t2pre*1000/dt;
                
                v1 = vars(t1:t2,1);
                v2 = vars(t1:t2,2);
                
                % animation stuff
%                 h = animatedline;
%                 for k = 1:35:length(v1)
%                     addpoints(h,v1(k),v2(k))
%                     drawnow
%                 end
%                 clearpoints(h)
                
                plot(v1,v2,'color',colors(c+1,:),'DisplayName',[mat2str(c) '0% ablated'])
                drawnow
                
                
            end
            
            legend('off')
            legend('show');
            set(get(gca,'YLabel'),'Rotation',0);
            
            saveas(f,['StateVariableDynamics/' obj.str '_indx' mat2str(figindex) '.png'])
            close all
            
        end
        
        function NetworkStatisticsDuringAblation(obj)
            
        end
        
    end
end

