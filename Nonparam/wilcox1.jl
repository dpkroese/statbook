using StatsBase, Distributions

x = [1.04, 3.30, 1.40, 1.53, 1.01, 1.02, 3.93, 1.61, 1.93, 2.25]
y = [1.00, 1.40, 1.30, 3.95, 0.08, 1.33, 0.66, 0.73, 1.49, 0.43]
m = length(x); n = length(y); N = m+n;
z = cat(x, y, dims=1)
ind = sortperm(z)
ranks = invperm(ind)
t = sum(ranks[1:10])

ET = m*(N+1)/2
varT = m*n*(N+1)/12
p = 2*(1 - cdf(Normal(0,1),(t - ET)/sqrt(varT)))

function pval(t)
   tot = 0
   a = 1:20
   for i = 0:2^20-1
      b = digits(i, base=2, pad=20)
      if sum(b)==10
        tot = tot + (sum(b .* a) >= t)
      end
   end
   return tot / binomial(20,10)
end

2*pval(t)

pooledV = ((m-1)var(x) + (n-1)var(y))/(N-2)
ttest = (mean(x) - mean(y))/sqrt(pooledV)/sqrt(1/m + 1/n) 
p1 = 2*(1 - cdf(TDist(N-2), ttest))