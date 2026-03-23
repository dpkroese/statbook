function sfran2_loglike(mu,eta_alpha,eta_gamma,eta,y,Xalpha,Xgamma)
   sigma2_alpha = exp.(eta_alpha)
   sigma2_gamma = exp.(eta_gamma)
   sigma2 = exp.(eta)
   n = length(y)
   Sigma = sigma2*diagm(ones(n)) .+ sigma2_alpha*(Xalpha*
   Xalpha') .+ sigma2_gamma*(Xgamma*Xgamma')
   l = -n/2*log(2*pi) - sum(log.(diag(cholesky(Sigma).L))) -
          0.5*(y .- mu)'*(Sigma\(y .- mu));
end
