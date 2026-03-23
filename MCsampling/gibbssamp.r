rm(list=ls());
cat("\f")
#dev.off();
f<-function(X1,X2){exp(-(X1*X2 + X1 +X2))*(X1 > 0 & X2 > 0)};
N = 10^4; #sample size
x = matrix(rep(0,2*N),,2);
x2=1;
for (i in 2:N){
  x1 = -log(runif(1))/(x2 +1);
  x2 = -log(runif(1))/(x1 +1);
  x[i,]=c(x1,x2);  
}
win.graph()
plot(x[,1],x[,2],pch=18,col="blue")