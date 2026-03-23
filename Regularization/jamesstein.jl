using StatsBase, LinearAlgebra,Plots
N = 100
#mu = rand(N)*100 .- 50; # not much difference
mu = rand(N)*0.1 .- 0.05 # a lot of difference
K = 10000
mse_js = zeros(K)
mse = zeros(K)
for k=1:K
   global mu,N
   X = mu + randn(N);
   mu_js = mean(X) .+ (1 - (N-3)/(var(X)*(N-2)))*(X .- mean(X))
   mse_js[k] = norm(mu_js - mu)^2
   mse[k] = norm(X - mu)^2
end
m1 = mean(mse_js)
m2 = mean(mse)
print("MSE JS =", m1,"   MSE MLE =  ",m2)
