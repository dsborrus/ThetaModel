function [meanperiod, meanAMP, brstbrstltratio] = psanalysis(rawspikes,dt,tmax)

tauavg = 1e2;
spikes = (rawspikes./tauavg)*1000;

u = mean(spikes);
o = std(spikes);

% spike threshold
spthre = u+2*o;
% burstlett threshold
bstlthre = 2;


% % % % % % % % % % % % % % % % % %

T = dt:dt:tmax;

L = length(spikes);

window = 300; % in frames
runwin = 0;

up = false;

spiketimes = [];
spikeV = [];



for z = 2:L
    runwin = runwin - 1;
    
    if runwin < 0
        
        if spikes(z) > bstlthre && spikes(z-1) < bstlthre
            up = true;
            runwin = window;
            vmax = spikes(z);
            vmaxT = z*dt;
        end
        
        
        if up
            
            if spikes(z)>vmax
                vmax = spikes(z);
                vmaxT = z*dt;
            end
            
            if spikes(z) < bstlthre && spikes(z-1) > bstlthre
                up = false;
                spikeV(length(spikeV)+1) = vmax;
                spiketimes(length(spiketimes)+1) = vmaxT; %#ok<*AGROW>
            end
            
        end
        
    end
    
end


meanperiod = mean(spiketimes);
meanAMP = mean(spikeV);

brstbrstltratio = length(find(spikeV>spthre))/length(spikeV);


figure; hold on;
plot(T,spikes);
plot(spiketimes,spikeV,'o','markersize',8);