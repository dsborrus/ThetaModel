%% Get synaptic depression nulcline
disp('Finding steady states for synaptic depression term')

tmax = 2e3;

Dresolution = 25;
Darray2 = linspace(Dlow,Dhigh,Dresolution);

ydrop = 0.1;
tauy = 5e3;
global MaxD %MaxD = 0.7;

spsarray = zeros(Dresolution,1);
yssarray = spsarray;

for i = 1:Dresolution
    [spsarray(i),yssarray(i)] = synctheta_v5_ynulcline(tmax,Darray2(i),ydrop,tauy); 
    disp([mat2str(i) ' of ' mat2str(Dresolution) ' done.'])
end

f2 = figure;
hold on;
plot(Darray/MaxD,Lowsps,'b'); hold on;
plot(Darray/MaxD,Highsps,'r');
plot(yssarray*MaxD,spsarray,'xg');
xlabel('y')
ylabel('Spikes per second')

save SDn_results