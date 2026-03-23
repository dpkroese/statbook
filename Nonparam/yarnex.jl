using StatsBase, Distributions

x = [15.0, 15.7, 14.8, 14.9, 13.0,  15.9, 
     17.0, 15.6, 15.8, 14.2, 16.2,  15.6,
     13.8, 17.6, 18.2, 15.0, 16.4,  15.0,
     15.5, 17.1, 16.0, 12.8, 14.8,  15.5]

y = [18.2, 17.2, 15.2, 15.6, 19.2, 16.2,
     16.8, 18.5, 15.9, 16.0, 18.0, 15.9,
     18.1, 15.0, 14.5, 15.2, 17.0, 14.9,
     17.0, 16.2, 14.2, 14.9, 16.9, 15.5]

# for bobbin 3:
 #    x = [14.8,15.8,18.2,16.0]
 #    y = [15.2,15.9,14.5,14.2]
m = length(x); n = length(y); N = m+n;
z = cat(x, y, dims=1)
ind = sortperm(z)
ranks = invperm(ind)
t = sum(ranks[1:m])

ET = m*(N+1)/2
varT = m*n*(N+1)/12
2*cdf(Normal(ET,sqrt(varT)),t)
#2*(1- cdf(Normal(ET,sqrt(varT)),t-0.5))

pooledV = ((m-1)var(x) + (n-1)var(y))/(N-2)
ttest = (mean(x) - mean(y))/sqrt(pooledV)/sqrt(1/m + 1/n) 
p1 = 2*cdf(TDist(N-2), ttest)