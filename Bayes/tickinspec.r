#tickinspec.r
rm(list=ls());
#dev.off();
source("kde.r");
source("gamrand.r");
c = 3.481048347*10^(19);
n = 10000;
x = 60;
p = cbind(1/3,1/3,1/3);
kk = matrix(0,1,n);
tt = matrix(0,1,n);
k = min(which(cumsum(p)> runif(1)));
for (i in 1:n){
  a = x+ 1;
  b = 2/k + 10*k; 
  t = gamrand(a,b);
  tt[i] = t;
  p = cbind(exp(-12*t),2^(x-1)*exp(-21*t), 3^(x-1)*exp(-92*t/3));
  p = p/sum(p);
  k = min(which(cumsum(p)> runif(1)));
  kk[i] = k;  
}
p1est = sum(kk == 1)/n  #estimate of posterior probability 1
p2est = sum(kk == 2)/n  #estimate of posterior probability 2
p3est = sum(kk == 3)/n  #estimate of posterior probability 3
p1est
p2est
p3est
result<-kde(tt,len=0,loop=1);
h = seq(0,4,0.01);
q = h^(60)*(3^(59)*exp(-92*h/3) + 2^59*exp(-21*h) + exp(-12*h))/c;
lines(h,q,lwd=2);
#plot(h,q,xlim=c(1,4),ylim=c(0,1.2),lwd=2,type='l',xlab=" t", ylab='f(t|x) = 60');