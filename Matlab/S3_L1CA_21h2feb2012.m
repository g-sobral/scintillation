%M = csvread('20h4Feb2012.txt');

S3_L1CA = M(M(:,2)==3 & M(:,3)==0,:);
freq = size(S3_L1CA,1)/(S3_L1CA(end,1)-S3_L1CA(1,1));

t0 = S3_L1CA(1,1);
time = S3_L1CA(:,1)-t0;

I = S3_L1CA(:,5);
Q = S3_L1CA(:,6);

A = zeros(size(time));
cIntensity = zeros(size(time));

for i=1:length(time)
  A(i) = sqrt(I(i)^2 + Q(i)^2);
  cIntensity(i) = A(i)^2;
end

% signal = Intensity(1:3000);
% n_signal = signal/mean(signal);
% cS4 = std(n_signal);

%plot(time(1:3000), 10*log10(signal))
Intensity = zeros(60,3000);
S4 = zeros(60,1);

for i=1:60
  Intensity(i,:) = cIntensity(1:3000);
  cIntensity = cIntensity(3000:end);
  
  signal = Intensity(i,:);
  n_signal = signal/mean(signal);
  S4(i) = std(n_signal);
end

signal1 = Intensity(34,:);
signal2 = Intensity(45,:);

n_signal1 = signal1/mean(signal1);
s4_signal1 = std(n_signal1)

n_signal2 = signal2/mean(signal2);
s4_signal2 = std(n_signal2)

figure
subplot(2,2,1)
plot(1/50:1/50:60, 10*log10(signal1))
axis([0 60 0 60])
xlabel('Time (s)')
ylabel('C/N_0 (dB-Hz)')
title(['Time series of intensity fluctuations S_{4}:' num2str(S4(34))])

subplot(2,2,2)
[psd_signal1,f1] = pwelch(n_signal1, [], [], [], 10);
semilogx(f1, 10*log10(psd_signal1))
axis([0.005 5 -50 20])


subplot(2,2,3)
plot(1/50:1/50:60, 10*log10(signal2))
axis([0 60 0 60])
xlabel('Time (s)')
ylabel('C/N_0 (dB-Hz)')
title(['Time series of intensity fluctuations S_{4}:' num2str(S4(45))])

subplot(2,2,4)
[psd_signal2, f2] = pwelch(n_signal2, [], [], [], 10);
semilogx(f2, 10*log10(psd_signal2))
axis([0.005 5 -50 20])
