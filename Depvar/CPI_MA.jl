using SparseArrays, LinearAlgebra, Optim, StatsBase,Plots,DelimitedFiles

function loglike_MA1(theta,y)
   # the log-likelihood function for MA(1)
   # input: theta = [psi sigma2]; y = data
   psi = theta[1]; sigma2 = theta[2];
   T = length(y);
   H = sparse(I,T,T) .+ psi*sparse(2:T,1:T-1,ones(T-1),T,T);
   HH = H*H';
   l = -T/2*log(2*pi*sigma2) - .5/sigma2*y'*(HH\y);
   return l
end

USCPI= readdlm("USCPI.csv")
y0 = USCPI[1]
y = USCPI[2:end]
Dely = y - [y0; y[1:end-1]]; # define the dependent variable
T = length(Dely)
theta0 = [0 ,var(Dely)]
f = theta -> -loglike_MA1(theta,Dely)
res = optimize(f,theta0)
thetahat = res.minimizer
psihat = thetahat[1];
l = loglike_MA1(thetahat,Dely);

Hhat = sparse(I,T,T) .+ psihat*sparse(2:T,1:T-1,ones(T-1),T,T);
uhat = Hhat\Dely;
