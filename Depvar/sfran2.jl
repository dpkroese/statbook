using Kronecker, LinearAlgebra, StatsBase, Optim

function sfran2_loglike(mu,eta_alpha,eta_gamma,eta,y,Xalpha,Xgamma)
sigma2_alpha = exp.(eta_alpha)
   sigma2_gamma = exp.(eta_gamma)
   sigma2 = exp.(eta)
   n = length(y)
   Sigma = sigma2*diagm(ones(n)) .+ sigma2_alpha*(Xalpha*Xalpha') .+
       sigma2_gamma*(Xgamma*Xgamma')
   l = -n/2*log(2*pi) - sum(log.(diag(cholesky(Sigma).L))) - 
          0.5*(y .- mu)'*(Sigma\(y .- mu));
end

y = [1.39, 1.29, 1.12, 1.16, 1.52, 1.62, 1.88, 1.87, 1.24, 1.18,
0.95, 0.96, 0.82, 0.92, 1.18, 1.20, 1.47, 1.41, 1.57, 1.65]
Xalpha = kronecker(diagm(ones(5)),ones(4,1))
Xgamma = kronecker(diagm(ones(10)),ones(2,1))
yhat = mean(reshape(y,4,5),dims=2)
theta0 = [mean(y), log(var(yhat)), log(var(y)/3), log(var(y)/3)];
f = theta -> -sfran2_loglike(theta[1],theta[2],theta[3], theta[4],y,Xalpha,Xgamma)
res = optimize(f,theta0);
thetahat = res.minimizer