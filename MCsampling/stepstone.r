rm(list=ls());
dev.off();
n = 101;
P=t(matrix(c(0,0.2,0,0.3,0.5,0,0.5,0,0.5,0,0,0,0.3,0,0,0.7,0,0,0,0,0,0,0,1,0,0,0,0.8,0,0.2,0,0,0.1,0,0.9,0),6,6));
x = matrix(0,1,n);
x[1]=1;
tot = matrix(0,1,6);
tot[1] = 1;
for (t in 1:(n-1))
{
  x[t+1] = min(which(cumsum(P[x[t],])>runif(1)))
  tot[x[t+1]] = tot[x[t+1]] + 1;
}
plot((0:(n-1)),x,pch=18,col="blue")
lines((0:(n-1)),x,col="blue")
tot/n
f = t(Null(t((diag(1,6)-t(P)))));
f = f/sum(f)