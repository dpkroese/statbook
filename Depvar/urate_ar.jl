urate = [5.7, 5.8, 5.7, 5.9, 5.9, 6.1, 6.1, 5.8, 5.7, 
5.6, 5.4, 5.4, 5.3, 5.1, 5.0, 5.0, 4.7, 4.6, 4.6, 4.4, 
4.5, 4.5, 4.7, 4.8, 5.0, 5.3, 6.0, 6.9, 8.3, 9.3, 9.6, 
9.9, 9.8, 9.6, 9.5, 9.6, 9.0, 9.0,  9.1, 8.7]
y = urate[3:end]
T = length(y)
X = [ones(T,1) urate[2:end-1]]
rhohat = (X'*X)\(X'*y)
yhat1 = X*rhohat               # fitted values
uhat = y-yhat1                 # residuals
sig2hat = uhat'*uhat/T

Z = [ones(T,1) urate[2:end-1] urate[1:end-2]]
phihat = (Z'*Z)\(Z'*y)
vhat = y-Z*phihat;

using Plots
t = 2002.5:.25:2011.75
plot(t,X*rhohat)
plot!(t,y)

plot(t,Z*phihat)
plot!(t,y)

using StatsBase
function acov(x,s)
    return sum((x[1:end-s] .- mean(x)).*(x[s+1:end] .- mean(x)))/
   (length(x) - s - 1);
end
acov(uhat,1)/acov(uhat,0)