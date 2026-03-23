rm(list=ls());
c = (2*pi)^(3/2);
#H <-function(x,y=X) {x<-c*sqrt(abs(matrix(rowSums(Z,1),nrow(Z),))); y} ;
H <-function(x,y=x) {x<-matrix(c*sqrt(abs(rowSums(x,1))),nrow(x),); y} ;

N = 10^6;
Z = matrix(c(rnorm(N),rnorm(N),rnorm(N)),ncol=3); X = H(Z);
mX = mean(X); sX=sd(as.numeric(X));
R=1.96*sX/sqrt(N);
output <- sprintf("a=%g, CI=(%g,%g)",mX,mX-R,mX+R)
print(output)