using Distributions, Plots
z = [40, 5, 1, -7, 15, 3, 12, -6, 2, 16]
ind = sortperm(abs.(z))
ranks = invperm(ind)
t = sum(ranks .* (z .> 0))

function pval(t)
   tot = 0
   a = 1:10
   for i = 0:2^10-1
      b = digits(i, base=2, pad=10)
      tot = tot + (sum(b .* a) >= t)
   end
   return tot / 2^10
end

pval(t)
tt = 0:0.1:55
plot(tt, pval.(tt))