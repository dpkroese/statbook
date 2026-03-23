using SparseArrays, LinearAlgebra, Optim, Plots,DelimitedFiles
function AR1_loglike(rho,y,X)
  T = length(y);
  H = sparse(I,T,T) .- rho*sparse(2:T,1:T-1,ones(T-1),T,T);
  HH = H'*H;
  betahat = (X'*HH*X)\(X'*HH*y);
  e = y-X*betahat;
  sigma2hat = e'*HH*e/T;
  l = -T/2*log(2*pi*sigma2hat) - .5/sigma2hat*e'*HH*e;
return l, betahat,sigma2hat
end

ads = readdlm("ads.csv",',')
y = ads[:,1]
T = length(y)
X = [ones(T,1) ads[:,2]]
f(rho)  = -(AR1_loglike(rho,y,X)[1])

res = optimize(f,0.1,1)
rhohat = res.minimizer
l, betahat, sigma2hat = AR1_loglike(rhohat,y,X)

e = y - X*betahat
H = sparse(I,T,T) .- rhohat*sparse(2:T,1:T-1,ones(T-1),T,T);
u = H*e;
scatter(u)
