using Distributions, Random, Statistics
mu = 3; sig = 0.5
# true parameters
alpha = 0.05; n = 10; mu_count = 0; sig_count = 0
for k in 1:100
   x = randn(n)*sig .+ mu  # draw the iid sample
   mu_est = mean(x)  # estimate mu
   sig_est = std(x)  # estimate sigma
   tq = quantile(TDist(n-1),1-alpha/2)
   mu_lo = mu_est - tq*sig_est/sqrt(n)  # low bound CI for mu
   mu_hi = mu_est + tq*sig_est/sqrt(n)  # upper bound
   cq1 = quantile(Chisq(n-1),1-alpha/2)
   cq2 = quantile(Chisq(n-1),alpha/2)
   sig_lo = (n-1)*sig_est^2/cq1; # lower bound CI for sigma
   sig_hi = (n-1)*sig_est^2/cq2;  # upper bound
   global mu_count = mu_count + (mu_lo < mu < mu_hi)
   global sig_count = sig_count + (sig_lo < sig^2 < sig_hi)
end
println(mu_count, "   ", sig_count) 