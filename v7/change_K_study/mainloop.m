% loop to run param study to range K

clear; close all; clc

t1=tic;

iu_min = -.00110;
iu_max = -.00080;

resolution = 40;

tmax = 7e5;
dt = 0.5;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

iurange = linspace(iu_min,iu_max,resolution);

mp = zeros(resolution,1);
mA = mp;
mbblr = mp;

for k = 1:resolution
    t2 = tic;
    spikes = synctheta_v7(tmax,iurange(k),4,0,0);
    
    [mp(k), mA(k), mbblr(k)] = psanalysis(spikes,dt,tmax);
    
    disp([mat2str(k) ' out of ' mat2str(resolution) ' done.'])
    toc(t1)
    disp(['Estimated time remaining is ' mat2str(round((toc(t2)/60)*(resolution-k))) ' minutes.'])
    
end

save AAA_data

figure;
subplot(3,1,1);
plot(iurange,mA); title('amplitude')
subplot(3,1,2);
plot(iurange,mp); title('period')
subplot(3,1,3);
plot(iurange,mbblr); title('burst:burstlett ratio')
xlabel('iu values');


toc
