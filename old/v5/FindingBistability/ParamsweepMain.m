%% Sweep D to find bistability
disp('Beginning parameter sweep of D to find bistability')
t1 = tic;

Dresolution = 59;
Dlow = 0.01;
Dhigh = 2.5;

Darray = linspace(Dlow,Dhigh,Dresolution);
Lowsps = zeros(Dresolution,1);
Highsps = Lowsps;

for i = 1:length(Darray)
    % synaptic depression off
    t2 = tic;
    
    % low state
    params.istate = 1;
    params.bumpit = 0;
    params.silence = 1;
    Lowsps(i) = synctheta_v5_PS(Darray(i),params);

    % high state
    params.istate = 2;
    params.bumpit = 1;
    params.silence = 0;    
    Highsps(i) = synctheta_v5_PS(Darray(i),params);
    
    if i~=Dresolution
        disp(['We have finished ' mat2str(i) ' out of ' mat2str(Dresolution) ...
            ' simulations.'])
        disp(['ETA is ' mat2str(round( toc(t2)*(Dresolution-i),2)) ' seconds.' ...
            ' or ' mat2str(round( toc(t2)*(Dresolution-i)/60)) ' minutes.'])
        
    end
end

f = figure;
subplot(3,3,[1 6])
plot(Darray,Lowsps,'b'); hold on;
plot(Darray,Highsps,'r');
xlabel('D')
ylabel('Spikes per second')

subplot(3,3,[7 9]); hold on;
str = ['isig1 = ' mat2str(params.isig1) '. iu1 = ' mat2str(params.iu1) ...
    '. isig2 = ' mat2str(params.isig2) '. iu2 = ' mat2str(params.iu2) ...
    '. n1 = ' mat2str(params.n1) '. n2 = ' mat2str(params.n2) '. prob = ' mat2str(params.prob) ...
    ];
annotation('textbox',[.2 .08 .1 .1],'String',str,'FitBoxToText','on');

save PSresults.mat

disp(['Total elapsed time for PS is ' mat2str(round(toc(t1),2)) ' seconds'])

