using Random, Statistics
xbar_obs = -0.7; s_obs = 0.4; t_obs = 2*xbar_obs/s_obs
N = 10^5;
count = 0;
for i in 1:N
    x = randn(4);
    xbar = mean(x); s = std(x); t = 2*xbar/s;
    global count = count + (t <= t_obs);
end
phat = count/N  # estimated p-value