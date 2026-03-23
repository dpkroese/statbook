rm(list=ls());
source("kde.r");
n = 100;
k = 5000;
set.seed(123456);
xorg =  tan(pi*(0.5 - matrix(runif(n),,n)));# original data
medxorg =  median(xorg);
meanxorg = mean(xorg);
x = matrix(0,1,n);
mx = matrix(0,1,k);
for (i in 1:k)
{
  ind = ceiling(n*matrix(runif(n),,n)); # draw random indices
  x = xorg[ind];# resampling the data (R)
  #x = tan(pi*(0.5 - matrix(runif(n),,n))); % sampling the data (S)
  mx[i] = median(x)
  #mx[i] = mean(x);
}
result<-kde(mx,2^7,len=3);
bandwidth<-result[[1]]
density<-result[[2]]
xmesh<-result[[3]]
win.graph()
plot(smooth.spline(xmesh,density),type='l',col='blue')