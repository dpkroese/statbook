# sfran2.r
source("sfran2_loglike.r");
y = matrix(c(1.39, 1.29, 1.12, 1.16, 1.52,
     1.62, 1.88, 1.87, 1.24, 1.18,
     0.95, 0.96, 0.82, 0.92, 1.18,
     1.20, 1.47, 1.41, 1.57, 1.65));
Xalpha=kronecker(diag(5),matrix(rep(1,4),,1));
Xgamma=kronecker(diag(10),matrix(rep(1,2),,1));
f<-function(theta){
  res<- sfran2_loglike(theta[1],theta[2],theta[3],theta[4],y,Xalpha,Xgamma);
  return (-res);
}
yhat = colMeans(matrix(y,4));
theta0 = c(mean(y), log(var(yhat)), log(var(y)/3),log(var(y)/3));
thetahat=optim(par=theta0, fn=f, gr = NULL, method = c("L-BFGS-B"),lower = -Inf, upper = Inf,control = list(), hessian = T)$par;
thetahat[1]
exp(thetahat[2:4])