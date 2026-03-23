rm(list=ls());
source("kde.r")
h = 0.1; 
h2 = h^2; 
c = 1/sqrt(2*pi)/h;
phi<-function(x,y){exp((-(x-y)^2/(2*h2)))}
f <-function(x) {exp(-x)*(x>=0);}
N = 10^4;
x = matrix(c(-log(runif(N))));
xx = matrix(seq(-0.5,6,0.01),1)
phis = matrix(0,1,length(xx))
for (i in 1:N)
{
  phis = phis + phi(xx,x[i]);
}
phis = c*phis/N;
win.graph();
plot(xx,phis,type='l',col="red");
output=list("bandwidth","density","xmesh")
len=length(output);
#[bandwidth,density,xmesh] = kde(x,2^12,0,max(x));
result<-kde(x,2^12,0,max(x),len)
bandwidth<-result[[1]]
density<-result[[2]]
xmesh<-result[[3]]
idx <- which(xmesh <= 6);
lines(xmesh[idx],density[idx],col="blue")
lines(xx,f(xx),col="black")
