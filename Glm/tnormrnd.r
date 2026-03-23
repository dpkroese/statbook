# tnormrnd.r
tnormrnd<-function(mu,sigma2,a,b){
  N = length(mu);
  sigma = sqrt(sigma2);
  u = matrix(runif(N));
  p1 = pnorm((a-mu)/sigma);
  p2 = pnorm((b-mu)/sigma);
  C = qnorm(p1+(p2-p1)*u);
  draw = mu +  sigma*C;
}
