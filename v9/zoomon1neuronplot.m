function zoomon1neuronplot(Parameters,Options,t,y,m,n,sij,spikes,raster,rr,Ihistory,sijhistory,theta)

nplots = 7;

%rr = 100;

figure('Position',[100 500 1500 1200]);

ax3 = subplot(nplots,1,1);
plot(t/1000,y(:,rr))
title(['Example neuron`s synaptic depression (y)'])

ax5 = subplot(nplots,1,2);
plot(t/1000,m(:,rr))
title(['Example random neuron`s fast synaptic depression term (m)'])

ax6 = subplot(nplots,1,3);
plot(t/1000,n(:,rr))
title(['Example neuron`s slow synaptic depression term (n)'])

ax4 = subplot(nplots,1,4);
plot(t/1000,sij(:,rr))
title(['Example neuron`s synaptic current output (s)'])

ax7 = subplot(nplots,1,5);
plot(t/1000,sijhistory)
title(['Example neuron`s synaptic current input'])

ax8 = subplot(nplots,1,6);
plot(t/1000,Ihistory)
xlabel('Time (s)')
title(['Example neuron`s total applied current (I+delta*syncurrent*syndepression)'])

ax9 = subplot(nplots,1,7);
plot(t/1000,theta(:,rr)); title('theta')

linkaxes([ax3 ax5 ax4 ax6 ax5 ax7 ax8 ax9],'x')

end