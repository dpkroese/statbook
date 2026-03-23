%kolsmirweib.m
clear all
rng(1234);
alpha = 1.5;
%generate data 
N = 100;
U = rand(1,N);
x = (-log(U)).^(1/alpha);
y = sort(1 - exp(-x));
i=1:N;
%plot empirical cdf
stairs([0,y],[0,i/N],'r'), hold on, line([0,1],[0,1]);
dn_up = max(abs(y-i/N));dn_down = max(abs(y-(i-1)/N));
dn = max(dn_up, dn_down);

%approximate p-value
j = -40:1:40;
p1 = 1 - sum((-1).^j.*exp(-2*(j.*sqrt(N)*dn).^2))

%or use matlab Statistics toolbox function kstest
[h,p2,ksstat,cv] = kstest(y,[y',y']);
p2

%or use Monte Carlo simulation
K = 1000;
for k=1:K
    i=1:N;
    y = sort(rand(1,N));
    DN(k) = max( max(abs(y-i/N)), max(abs(y-(i-1)/N)));
end
p  = sum(DN >= dn)/K

