using Random, Plots, StatsBase , StatsPlots
Random.seed!(1234); 
alpha = 1.5;
N = 100;
U = rand(N);
x = (-log.(U)).^(1/alpha); # generate data 
y = sort(1 .- exp.(-x));
i = 1:N;

ecdfplot(y,legend=false)  # empirical cdf
plot!([0,1],[0,1])
dn_up = maximum(abs.(y -i/N));
dn_down = maximum(abs.(y -(i .- 1)/N));
dn = max(dn_up, dn_down);

# Use MC simulation to obtain the p-value
K = 10000;
DN = zeros(K);
for k in 1:K
    local i = 1:N;  # a global i already exists
    local y = sort(rand(N)); # same for y
    DN[k] = max( maximum(abs.(y-i/N)), maximum(abs.(y-(i .- 1)/N)));
end
p  = sum(DN .>= dn)/K

j = -40:1:40;
p1 = 1 - sum((-1).^j.*exp.(-2*(j.*sqrt(N)*dn).^2))

#Q7.6a