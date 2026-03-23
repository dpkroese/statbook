N= 1000;
rho = 0.8;
Sigma<-rbind(c(1,rho), c(rho,1));
B=t(chol(Sigma))
x=B%*%matrix(rnorm(N),2)
plot(x[1,],x[2,],pch=18,col="blue")