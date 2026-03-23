include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Plots, Distributions, .ThetaKDE
n = 10000;
x = 60; #number of tickets
p = [1/3, 1/3, 1/3] #initial value
kk = zeros(Int64,n);
tt = zeros(n);
k = minimum(findall(cumsum(p) .> rand()));
for i in 1:n
    a = x+ 1;
    b = 2/k + 10*k;
    t = rand(Gamma(a,1/b));
    tt[i] = t;
    global p = [exp(-12*t), 2.0^(x-1)*exp(-21*t),  3.0^(x-1)*exp(-92*t/3)];
    p = p/sum(p);
    global k = minimum(findall(cumsum(p) .> rand()));
    kk[i] = k;
end
p1est = sum(kk .== 1)/n #estimate of post. prob. 1
p2est = sum(kk .== 2)/n #estimate of post. prob. 2
p3est = sum(kk .== 3)/n #estimate of post. prob 3
xmesh, density, bw = kde(tt)
plot(xmesh,density)
f1(t) = t^60
f2(t) = 3.0^59*exp(-92*t/3)
f3(t) = 2.0^59*exp(-21*t)
f4(t) = exp(-12*t)
tickf(t) = f1(t)*(f2(t) + f3(t)+f4(t))/3.481048347e19
plot!(xmesh,tickf.(xmesh))