include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Random, Plots, QuadGK, .ThetaKDE
N = 10^5;   # sample size
f(x) = x^2*exp(-x^2 + sin(x)); # unnormalized target pdf
g(x) = exp(-abs(x))/2;         # proposal pdf
alpha(x,y) = min(f(y)*g(x)/(f(x)*g(y)), 1); # accept. prob.
x = 0; xx = zeros(N);
for t in 2:N
    global x
    y = -log(rand())*(2*(rand() < 1/2) - 1);  # proposal
    rand() < alpha(x,y) ? x = y : nothing
    xx[t] = x;
end
jx = xx[1:N] + randn(N)*0.05;
xmesh,density,bw = kde(jx);
plot(xmesh,density)  # plot the kde of the data
c = quadgk(f,-5,5)[1]; # determine the normalization constant
tt = -4:0.1:4; 
plot!(tt,f.(tt)/c) # plot the target pdf 