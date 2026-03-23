using SparseArrays, LinearAlgebra, Optim, StatsBase,Plots,DelimitedFiles
function loglike_ARMA11(psi,X,y)
   T = length(y)
   H = sparse(I,T,T) + psi*sparse(2:T,1:T-1,ones(T-1),T,T)
   HH = H*H'
   rhohat = (X'*(HH\X))\(X'*(HH\y))
   uhat = y-X*rhohat
   sigma2hat = uhat'*(HH\uhat)/T
   l = -T/2*log(2*pi*sigma2hat) - .5/sigma2hat*uhat'*(HH\uhat);
   return l, rhohat, sigma2hat
end

USCPI= readdlm("USCPI.csv")
y0 = USCPI[1]
y = USCPI[2:end]
T = length(y)
X = [ones(T,1) [y0; y[1:end-1]]]
f = psi ->  -loglike_ARMA11(psi,X,y)[1]
res = optimize(f,-1,1)
psihat = res.minimizer
l, rhohat, sigma2hat = loglike_ARMA11(psihat,X,y)
