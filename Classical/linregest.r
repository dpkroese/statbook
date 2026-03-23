rm(list=ls())
N=100;
x=matrix((1:N)/N,N,1)
beta=matrix(c(6,13),2);
sigma=2;
A=matrix(c(matrix(1,100),x),,2)
y=A%*%beta + sigma*matrix(c(rnorm(N)),,1);
plot(x,y,pch=18,col="blue")
betahat=solve(t(A)%*%(A))%*%(t(A)%*%y)
sigmahat=norm(y-A%*%betahat,"f")/sqrt(N)
