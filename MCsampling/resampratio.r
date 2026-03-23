rm(list=ls());
#dev.off();
source("kde.r");
n = 100; 
k = 50000;
est = matrix(0,1,k);
xorg = 11 + 5*matrix(rnorm(n),,n);
yorg = matrix(runif(n),,n)*xorg;
x = matrix(0,1,n);
y = matrix(0,1,n);
estorg =  mean(xorg)/mean(yorg);
est = matrix(0,1,k);
for (i in 1:k)
{
  ind = ceiling(n*matrix(runif(n),,n));
  x = xorg[ind];
  y = yorg[ind];
  est[i] = mean(x)/mean(y);
}
result<-kde(est,len=0);