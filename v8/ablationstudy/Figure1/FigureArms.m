classdef FigureArms < handle
    %FIGUREPARTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Parameters
        Options
        str
        A
        Aexists
        
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
            obj.Parameters.tmax = 100000;
            obj.Options.Doablate    = 0;

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
                

            % plot
            tstart = 1500;
            tend = obj.Parameters.tmax;
            dt = obj.Parameters.dt;
            tauavg = obj.Parameters.tauavg;

            tw = tstart/dt:tend/dt;

            f = figure; hold on;
            plot(out.t(tw)/1000,(out.spikes(tw)./tauavg)*1000)
            xlim([out.t(tw(1))/1000 out.t(tw(end))/1000])
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
            obj.Parameters.tmax = 2250000; %(in ms)
            
            % overwriting ablation params as well!!!
            obj.Parameters.t2ablat  = 250;     % When to ablate neurons (every XXX seconds)
            obj.Parameters.N2k_w1p  = 5;       % number of neurons to kill with one pulse
            
            % remember to ablate
            obj.Options.Doablate = 1;

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
                

            % plot
            tstart = 1000;
            tend = obj.Parameters.tmax;
            dt = obj.Parameters.dt;
            tauavg = obj.Parameters.tauavg;

            tw = tstart/dt:tend/dt;

            f = figure; hold on;
            plot(out.t(tw)/1000,(out.spikes(tw)./tauavg)*1000)
            xlim([out.t(tw(1))/1000 out.t(tw(end))/1000])
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
            plot(obj.outablate.t(tw)/1000,(obj.outablate.spikes(tw)./tauavg)*1000,'linewidth',2)
            xlim([obj.outablate.t(tw(1))/1000 obj.outablate.t(tw(end))/1000])
            %ylim([0 20])
            %plot(obj.Parameters.t2ablat:obj.Parameters.t2ablat:obj.Parameters.tmax/1000,18,'r.','markersize',30)
            %title('Network Activity (intact)')
            xlabel('Time (s)')
            ylabel('Spikes/neuron/second')
            set(gca,'fontsize',23)
            
            % hardcoding in what folder to save to
            switch chunknumber
                case 2
                    ylim([0 15])
                    saveas(f,['10percent_ablated/' obj.str '.png'])
                case 4
                    ylim([0 15])
                    saveas(f,['20percent_ablated/' obj.str '.png'])
                case 6
                    ylim([0 15])
                    saveas(f,['30percent_ablated/' obj.str '.png'])
                case 8
                    ylim([0 15])
                    saveas(f,['40percent_ablated/' obj.str '.png'])
            end
                    
            close all
            
            
        end
        
    end
end

