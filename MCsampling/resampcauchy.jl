include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Random, Plots, StatsBase, .ThetaKDE

N = 100; K = 5000
#Random.seed!(123)
xorg = tan.(pi*(0.5 .- rand(N))) # original data
medxorg = median(xorg); meanxorg = mean(xorg)
x = zeros(N); mx = zeros(K)
for i in 1:K
    ind = ceil.(Int64,N*rand(N)) # draw random indices
    x = xorg[ind]; # resampling the data (R)
    # x = tan(pi*(0.5 - rand(1,N))) # sampling new data (S)
    mx[i] = median(x)
    # mx(i) = mean(x);
end

xmesh,density,bw = kde(mx,res=true)
plot(xmesh,density)

x = zeros(N); mx = zeros(K)
for i in 1:K
    ind = ceil.(Int64,N*rand(N)) # draw random indices
    x = xorg[ind];
    # resampling the data (R)
    # x = tan(pi*(0.5 - rand(1,N))) # sampling new data (S)
    #mx[i] = median(x)
     mx[i] = mean(x);
end

xmesh, density, bw=kde(mx,res=true)
plot!(xmesh,density)