%M = csvread('20h4Feb2012.txt');

satellites = unique(M(:,2));
freq = zeros(size(satellites));

for i=1:length(satellites)
  L1CA = M(M(:,2)==satellites(i) & M(:,3)==0);
  if(size(L1CA,1)>0)
    freq(i) = size(L1CA,1)/(L1CA(end,1)-L1CA(1,1));
  else
    freq(i) = 0;
  end
end

[satellites freq]