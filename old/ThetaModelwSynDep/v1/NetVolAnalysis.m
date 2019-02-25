function [gradeMa,gradeMu] = NetVolAnalysis(NetworkActivity,dt)
%NETVOLANALYSIS Grades the network activity on how accurate it is to
%desired shape
% uses fft right now

% NetworkActivity = spikes

L = length(NetworkActivity);     % Signal length
Fs = 1/dt;                       % Sampling Frequency

% Compute the Fourier transform of the signal
Y = fft(NetworkActivity);      

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum
% P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f and plot the single-sided amplitude
% spectrum P1.
f = Fs*(0:(L/2))/L;

% Plot
if 0
figure
plot(f,P1)
title('Signle-Sided Amplitude Spectrum of Network Actiovity')
xlabel('f (Hz)')
ylabel('P1(f)')
xlim([0 1])
end

% Grade system
% grading for max power between frequencies .2(1/5 seconds) and 
% 0.083 (1/12 seconds)
minf = 0.083; maxf = 0.2;
f_indx = f>minf & f<maxf;
fimp = P1(f_indx); % important frequencies

gradeMa = max(fimp);
gradeMu = mean(fimp);

end

