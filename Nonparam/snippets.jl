using Distributions
1 - cdf(Binomial(10),7)
1 - cdf(TDist(9),27/sqrt(209))
n = 10
t =  10+4 + 1 + 8 + 3 + 7 + 2 + 9
et = n*(n+1)/4
sdt = sqrt(n*(n+1)*(2*n+1)/24)
pval = 1 - cdf(Normal(0,1),(t - et)/sdt)
print(pval)
using Distributions
n = 10; t = 44
et = n*(n+1)/4
sdt = sqrt(n*(n+1)*(2*n+1)/24)
pval = 1 - cdf(Normal(et,sdt), t)
print(pval)