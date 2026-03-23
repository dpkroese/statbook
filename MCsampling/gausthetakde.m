%gausthetakde.m
clear all
h = 0.1; h2 = h^2; c = 1/sqrt(2*pi)/h;
phi = @(x,x0) exp(-(x-x0).^2/(2*h2)); % unscaled Gaussian kernel 
f = @(x) exp(-x).*(x >= 0);  % True pdf
N = 10^4;   % sample size
x = -log(rand(N,1)); % generate the data
xx = [-0.5:0.01:6];  % plot range
phis = zeros(1,numel(xx));
for i = 1:N
    phis = phis + phi(xx,x(i));
end
phis = c*phis/N;
hold on
plot(xx,phis,'r');        % plot Gaussian KDE
[bandwidth,density,xmesh] = kde(x,2^12,0,max(x));
idx = find(xmesh <= 6);  
plot(xmesh(idx),density(idx),'b')  % plot theta KDE
plot(xx,f(xx),'k');      % plot true pdf
