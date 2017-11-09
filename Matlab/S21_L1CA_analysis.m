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

max(S4)