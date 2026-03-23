using Distributions
z =   [40 , 5,  1 , -7, 15 , 3 , 12 , -6 , 2 , 16]
ind = sortperm(abs.(z))
ranks = invperm(ind)
t = sum(ranks .* (z .> 0))

function pval(a,K,t)
   tot = 0
   for i=1:K
     b = rand(10) .< 0.5
     tot = tot + (sum(b .* a) >= t)
   end
return tot/K
end

K = 1e6
a = 1:10
pval(a,K,t)

function pvale(a,t)
tot = 0
n = 10
for i = 0:2^n-1
   b = digits(i, base=2, pad=n)
   tot = tot + (sum(b .* a) >= t)
end
return tot/2^n
end

pvale(a,t)