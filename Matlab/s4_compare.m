close all; clear all; clc;

load RXModelAnalysis.mat

s4 = zeros(1,length(Intensity(1,:)));

for i=1:length(Intensity(1,:))
  signal = Intensity(:,i);
  n_signal = signal/mean(signal);
  s4(i) = std(n_signal);
end

sum(s4-sS4s)
sum(s4-eS4)
sum(s4-s_four)