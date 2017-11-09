close all; clear all; clc;

load RXModelAnalysis.mat

signal1 = Intensity(:,1221);
signal2 = Intensity(:,1256);

n_signal1 = signal1/mean(signal1);
s4_signal1 = std(n_signal1)

n_signal2 = signal2/mean(signal2);
s4_signal2 = std(n_signal2)

figure
subplot(2,2,1)
plot(0.1:0.1:300, 10*log10(signal1))
axis([0 300 0 70])
xlabel('Time (s)')
ylabel('C/N_0 (dB-Hz)')
title(['Time series of intensity fluctuations S_{4}:' num2str(eS4(1221))])

subplot(2,2,2)
[psd_signal1,f1] = pwelch(n_signal1, [], [], [], 10);
semilogx(f1, 10*log10(psd_signal1))
axis([0.005 5 -50 20])


subplot(2,2,3)
plot(0.1:0.1:300, 10*log10(signal2))
axis([0 300 0 70])
xlabel('Time (s)')
ylabel('C/N_0 (dB-Hz)')
title(['Time series of intensity fluctuations S_{4}:' num2str(eS4(1256))])

subplot(2,2,4)
[psd_signal2, f2] = pwelch(n_signal2, [], [], [], 10);
semilogx(f2, 10*log10(psd_signal2))
axis([0.005 5 -50 20])
