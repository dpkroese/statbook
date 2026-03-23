rm(list=ls());
win.graph();
N=100; x=matrix((1:N)/N,N,1)
beta=matrix(c(6,13),2);
sigma=2;
X=matrix(c(matrix(1,N,1),x),ncol=2);
y=X%*%beta + sigma*matrix(rnorm(N),,1);
plot(x,y,pch=18,col="blue");
betahat=solve(t(X)%*%(X))%*%(t(X)%*%y);
sigmahat=norm(y-X%*%betahat,"f")/sqrt(N);
tquant=qt(0.975,N-2);
ucl = matrix(0,,N); lcl = matrix(0,,N);
rl = matrix(0,,N); el = matrix(0,,N); u=0;
for (i in 1:N)
{
  u = u + 1/N;
  a = matrix(c(1,u),2);
  rl[i] = t(a)%*%beta;
  ucl[i] = t(a)%*%betahat + tquant*norm(y-X%*%betahat,"f")*sqrt(t(a)%*%solve(t(X)%*%X)%*%a)/sqrt(N-2);
  lcl[i] = t(a)%*%betahat - tquant*norm(y-X%*%betahat,"f")*sqrt(t(a)%*%solve(t(X)%*%X)%*%a)/sqrt(N-2);
}
lines(x,rl,col='blue');
lines(x,ucl,col='blue');
lines(x,lcl,col='blue');
