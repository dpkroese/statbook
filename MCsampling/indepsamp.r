rm(list=ls());
source("kde.r");
#dev.off();
N = 10^5;
xx= matrix(0,1,N);
lam = 1;
f<-function(X){X<-X^2*exp(-X^2 + sin(X));X};
g<-function(X){X<-lam*exp(-abs(X)*lam)/2;X};
alpha<-function(X,Y){min(f(y)*g(x)/(f(x)*g(y)), 1)};
x = 0; xx[1] = x;
for (t in 2:(N+1))
{
  y = -log(runif(1))*(2*(runif(1) < 1/2) - 1);
  if (runif(1) < alpha(x,y))
    {
    x = y;
    }
xx[t] = x;
}
result<-kde(rbind(xx[1:(N+1)]),len=0);
c = integrate(f,-5,5);
tt=matrix(seq(-4,4,0.1),1);
lines(tt,f(tt)/as.numeric(c[1]),col='red')