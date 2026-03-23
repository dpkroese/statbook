%resampratio.m
clear all, clf
n = 100; %size of data
k = 50000; %resample size
est = zeros(1,k);
xorg = 11 + 5*randn(1,n); %original x data
yorg = rand(1,n).*xorg;   %original y data
x = zeros(1,n); y = zeros(1,n);
estorg =  mean(xorg)/mean(yorg);
est = zeros(1,k);
for i=1:k
    ind = ceil(n*rand(1,n)); % draw random indices
    x = xorg(ind); y = yorg(ind); % resampled  data
    est(i) = mean(x)/mean(y);
end
kde(est);
