include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Random, Plots, StatsBase, Distributions, .ThetaKDE
Random.seed!(123)
N = 100 #size of data
K = 5000 #resample size
est = zeros(K)
xorg = 11 .+ 5*randn(N); yorg = rand(N).*xorg; #orig. data
estorg = mean(xorg)/mean(yorg)
x = zeros(N); y = zeros(N);
est = zeros(K);
for i in 1:K
    ind = ceil.(Int64,N*rand(N)) # draw random indices
    local x = xorg[ind]; local y = yorg[ind]; # resampled data
    est[i] = mean(x)/mean(y);
end
xmesh,density,bandwidth = kde(est,res=true)
plot(xmesh,density)
cv = cov(hcat(xorg,yorg))
sigma2 = estorg^2*(var(xorg)/mean(xorg)^2 + var(yorg)/mean(yorg)^2 - 2*cv[1,2]/mean(xorg)/mean(yorg));
t = estorg-4*sqrt(sigma2/N):0.01: estorg+4*sqrt(sigma2/N);
z = pdf.(Normal(estorg,sqrt(sigma2/N)),t);
plot!(t,z) 