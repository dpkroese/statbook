%pvalsim.m
clear all
xbar_obs = -0.7; s_obs = 0.4;
t_obs = 2*xbar_obs/s_obs;
N = 10^5;
tot = 0;
for i=1:N
    x = randn(4,1);
    s = std(x); xbar = mean(x);
    t = 2*xbar/s;
    tot = tot + (t <= -3.5); 
end
phat = tot/N
cdf('t',-3.5,3)