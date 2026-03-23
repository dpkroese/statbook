using Plots, LinearAlgebra, StatsBase
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,3.8,4.8,1.7,-0.4,4.5,1.3,0.4,2.6,4,2.9,1.6]
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,8.05,3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04]
n = length(x)
xi = quantile(x,[.25,.5,.75])  # knots
K = length(xi)

# cubic polynomial and cubic spline
X = hcat(ones(n),x,x.^2,x.^3)
beta_cp = (X'*X)\(X'*y)
for k=1:K
    global X
    tmpX = (x .- xi[k]).^3;
    X =hcat(X,tmpX.*(tmpX.>0))
end
beta_cs = (X'*X)\(X'*y)

# natural cubic spline
R = [0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 1 1; 0 0 0 0 xi']
r = zeros(4)
A = vcat(hcat(X'*X,R'), hcat(R,zeros(4,4)))
beta_ns = A\vcat(X'*y,r)

xtilde = minimum(x):0.01:maximum(x)
ngrid = length(xtilde)
Xtilde = hcat(ones(ngrid),xtilde,xtilde.^2,xtilde.^3)
cp = Xtilde*beta_cp
for k=1:K
    global Xtilde
    tmp = (xtilde .- xi[k]).^3
    Xtilde = hcat(Xtilde,tmp.*(tmp.>0))
end
cspline = Xtilde*beta_cs
nspline = Xtilde*beta_ns[1:K+4]

plot(xtilde,cp,lw=2,color=:black,ls=:dash)
plot!(xtilde,cspline,lw=2,color=:blue,legend=false)
plot!(xtilde,nspline,lw=2,color=:red,legend=false)
scatter!(x,y,color=:black)
