%% Get synaptic depression nulcline
disp('Finding steady states for synaptic depression term')

params.istate = 3;
params.silence = 0;
params.bumpit = 0;

Dresolution = 50;
Darray2 = linspace(Dlow,Dhigh,Dresolution);

spsarray = zeros(Dresolution,1);
yssarray = spsarray;
yssn2_array = spsarray;

parfor i = 1:Dresolution
    %t1 = tic;
    [spsarray(i),yssarray(i),yssn2_array(i)] = synctheta_v6_ynulcline(params,Darray2(i)); 
    disp([mat2str(i) ' of ' mat2str(Dresolution) ' done.'])
    %disp(['Estimated time remaining = ' mat2str(round(toc(t1)*(Dresolution-1),2)) ' seconds. Or '...
    %    mat2str(round(toc(t1)*(Dresolution-1)/60)) ' minutes.'])
end

f2 = figure; subplot(3,3,[1 6])
hold on;
plot(Darray/params.MaxD,Lowsps,'b'); hold on;
plot(Darray/params.MaxD,Highsps,'r');
plot(yssarray,spsarray,'--','color',[0,0,1]);
plot(yssn2_array,spsarray,'--','color',[0,0.5,0]);
xlabel('y')
ylabel('Spikes per second')

subplot(3,3,[7 9]); hold on;
str = ['isig1 = ' mat2str(params.isig1) '. iu1 = ' mat2str(params.iu1) ...
    '. isig2 = ' mat2str(params.isig2) '. iu2 = ' mat2str(params.iu2) ...
    '. n1 = ' mat2str(params.n1) '. n2 = ' mat2str(params.n2) '. prob = ' mat2str(params.prob) ...
    '. ydrop = ' mat2str(params.ydrop) '. tauy = ' mat2str(params.tauy) ...
    ];
annotation('textbox',[.2 .08 .1 .1],'String',str,'FitBoxToText','on');

save SDn_results