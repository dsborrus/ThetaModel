%% Sweep D to find bistability
disp('Beginning parameter sweep of D to find bistability')
t1 = tic;

Dresolution = 10;
Dlow = .5;
Dhigh = 1;
tmax = 1e4;

Darray = linspace(Dlow,Dhigh,Dresolution);
Lowsps = zeros(Dresolution,1);
Highsps = Lowsps;

for i = 1:length(Darray)
    % synaptic depression off
    t2 = tic;
    
    % low state
    istate = 1;
    bumpit = 0;
    silence = 1;
    Lowsps(i) = synctheta_v5_PS(tmax,istate,bumpit,silence,Darray(i));

    % high state
    istate = 2;
    bumpit = 1;
    silence = 0;    
    Highsps(i) = synctheta_v5_PS(tmax,istate,bumpit,silence,Darray(i));
    
    if i~=Dresolution
        disp(['We have finished ' mat2str(i) ' out of ' mat2str(Dresolution) ...
            ' simulations.'])
        disp(['ETA is ' mat2str(round( toc(t2)*(Dresolution-i),2)) ' seconds.' ...
            ' or ' mat2str(round( toc(t2)*(Dresolution-i)/60)) ' minutes.'])
        
    end
end

f = figure;
plot(Darray,Lowsps,'b'); hold on;
plot(Darray,Highsps,'r');
xlabel('D')
ylabel('Spikes per second')

save PSresults.mat

disp(['Total elapsed time for PS is ' mat2str(round(toc(t1),2)) ' seconds'])

