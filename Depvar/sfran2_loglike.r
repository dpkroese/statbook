# sfran2_loglike.r
sfran2_loglike<-function(mu,eta_alpha,eta_gamma,eta,y,Xalpha,Xgamma){
  sigma2_alpha = exp(eta_alpha);
  sigma2_gamma = exp(eta_gamma);
  sigma2 = exp(eta);
  n = length(y);
  Sigma = sigma2*diag(n) + sigma2_alpha*(Xalpha %*% t(Xalpha))+sigma2_gamma*(Xgamma %*% t(Xgamma));
  l = -n/2*log(2*pi) - sum(log(diag(chol(Sigma)))) - .5*t(y-mu) %*% solve(Sigma,(y-mu));
  return (l);
}
