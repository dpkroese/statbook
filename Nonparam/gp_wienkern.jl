using LinearAlgebra, Plots
xx = 0:0.001:1
kappa(x, u) = min(x, u)
# sigma = 0.2;
# kappa(x,u) = exp(-(x-u)^2/(2*sigma^2))
n = length(xx)
K = zeros(n, n)
for i = 1:n
   for j = 1:n
      K[i, j] = kappa(xx[i], xx[j])
   end
end
U, S, V = svd(K);
B = U * diagm(sqrt.(S));
g = plot()
for i = 1:5
   global g
   yy = B * randn(n)
   g = plot!(xx, yy, legend=false)
   display(g)
end