%resampcauchy.m
clear all 
n = 100;
k = 5000;
rand('state',123456); %comment out when randomising the experiment
xorg =  tan(pi*(0.5 - rand(1,n))); % original data
medxorg =  median(xorg);
meanxorg = mean(xorg);
x = zeros(1,n); mx = zeros(1,k); 
for i=1:k
    ind = ceil(n*rand(1,n)); % draw random indices
    x = xorg(ind);  % resampling the data (R)
    % x = tan(pi*(0.5 - rand(1,n))); % sampling the data (S)
    mx(i) = median(x);
    % mx(i) = mean(x);
end
[bandwidth,density,xmesh]=kde(mx,2^7);
plot(xmesh,density)
