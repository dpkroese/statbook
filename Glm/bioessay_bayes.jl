include("bioassay.jl")
using Distributions
B = cholesky(V).L;
burnin = 500;
nloop = 10000 + burnin;
store_beta = zeros(nloop, 2);
nu = 5;   # df for the proposal
# log posterior density
logf(b) = (sum((y .- 1) .* (X * b) - log.((1 .+ exp.(-X * b)))))[1];
# log density of the t proposal
logprop(b) = (-0.5*(nu+2)*log(1 .+ 
             (b - betat)' * (V \ (b - betat)) / nu))[1];
beta = betat;     # initialize the chain
for i = 1:nloop
    global beta
    # candidate draw from the t proposal
    betac = betat + B * randn(2, 1) * sqrt(nu / rand(Gamma(nu / 2, 2)))
    rho = logf(betac) - logf(beta) + logprop(beta) - logprop(betac)
    exp(rho) > rand()  ? beta = betac : nothing
    store_beta[i, :] = beta'
end
store_beta = store_beta[burnin+1:end, :];  # discard the burnin
cov(store_beta)
mean(store_beta, dims=1)