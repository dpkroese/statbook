using  Distributions, Plots
N = 5000
z = randn(N) # test statistic
pi0 = 0.05
mu1 = 4
istore = zeros(Bool, N);
z = randn(N)
for i = 1:N
   if rand() < pi0
      z[i] = mu1 + 2 * randn()
      istore[i] = true
   end
end
pvals = 2 * (1 .- cdf.(Normal(0, 1), abs.(z)))
spvals = sort(pvals)  # sorted p-values
idx = sortperm(pvals) # the corresponding indices
q = 0.05
n = 300 # plot only first n sorted p-values
p1 = plot(spvals[1:n])
p1 = plot!(1:n, (1:n) * q / N)
maximum(findall(spvals .<= (1:N) * q / N))
padj = zeros(N)
c = 1
for i = N:-1:1
   global c
   padj[i] = min(c, spvals[i] * N / i)
   c = padj[i]
end
findall(padj .<= q)
p2 = plot(padj[1:n])
p2 = plot!(1:n, q * ones(n))
plot(p1, p2, layout=(2, 1))
