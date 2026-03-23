% tnormrnd.m
function [ draw ] = tnormrnd(mu,sigma2,a,b)
N = length( mu );
sigma = sqrt(sigma2);
u = rand(N,1);
p1 = normcdf((a-mu)./sigma);
p2 = normcdf((b-mu)./sigma);
C = norminv(p1+(p2-p1).*u);
draw = mu + sigma.*C;