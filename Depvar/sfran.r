# sfran.r
source("sfran_loglike.r");
y = matrix(c(22.6, 20.5, 20.8, 22.6, 21.2, 20.5,
             17.3, 16.2, 16.6, 21.4, 23.7, 23.2, 
             20.9, 22.2, 22.6, 14.5, 10.5, 12.3, 
             20.8, 19.1, 21.3, 17.4, 18.6, 18.6, 
             25.1, 24.8, 24.9, 14.9, 16.3, 16.6),ncol=3,byrow=TRUE);
f<-function(theta){
  res<- sfran_loglike(theta[1],theta[2],theta[3],y);
  return (-res);
}
yhat = matrix(rowMeans(y,2));
theta0 = c(mean(yhat), var(yhat), mean(apply(y,1,var)));
thetahat=optim(par=theta0, fn=f, gr = NULL, method = c("L-BFGS-B"),lower = -Inf, upper = Inf,control = list(), hessian = T)$par;
print(thetahat);