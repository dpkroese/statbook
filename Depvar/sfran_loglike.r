# sfran_loglike.r
sfran_loglike<-function(mu,sigma2_mu,sigma2,y){
  d = dim(y)[1];
  ni = dim(y)[2];
  Sigmai = sigma2*diag(ni) + sigma2_mu*matrix(rep(1,ni*ni),nrow=ni);
  l = -(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai));
  for (i in 1:d){
  yi = matrix(y[i,]);
  l = l - .5*t(yi-mu) %*% solve(Sigmai,(yi-mu));    
  }
  return (l);
}