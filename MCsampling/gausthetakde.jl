include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Random, Plots, .ThetaKDE
h = 0.1; h2 = h^2; c = 1/sqrt(2*pi)/h; # constants

phi(x,x0) = exp(-(x -x0)^2/(2*h2)) # unscaled kernel
f(x) = x >= 0 ? exp(-x) : 0 # true pdf  

N = 10^4 # sample size
x = -log.(rand(N)) # generate the data

# Determine theta KDE
xmesh, density, bw = kde(x);
# Determine Gaussian KDE
phis = zeros(length(xmesh))
for i=1:N
   global phis = phis + phi.(xmesh,x[i])
end
phis = c*phis/N
plot(xmesh,phis)       # Gaussian KDE
plot!(xmesh,density)   # theta KDE
plot!(xmesh,f.(xmesh)) # true pdf