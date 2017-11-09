% close all; clear all; clc;
% M = csvread('21h2Feb2012.txt');

satellites = unique(M(:,2));
freq = zeros(size(satellites));
S4 = zeros(size(satellites));

for i=1:length(satellites)
  L1CA = M(M(:,2)==satellites(i) & M(:,3)==0,:);
  
  if(size(L1CA,1)>0)
    freq(i) = size(L1CA,1)/(L1CA(end,1)-L1CA(1,1));
    
    I = L1CA(:,5);
    Q = L1CA(:,6);

    A = zeros(size(L1CA,1),1);
    Intensity = zeros(size(L1CA,1),1);

    for j=1:size(L1CA,1)
      A(j) = sqrt(I(j)^2 + Q(j)^2);
      Intensity(j) = A(j)^2;
    end
    
    n_Intensity = Intensity/mean(Intensity);
    S4(i) = std(n_Intensity);
    
  else
    freq(i) = 0;
  end   
end

sortrows([satellites freq S4],3)