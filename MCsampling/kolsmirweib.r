rm(list=ls());
set.seed(123)
alpha = 1.5;
#generate data
N = 100;
U = matrix(runif(N),,N)
x = matrix(c(as.vector(-log(U))^(1/alpha)),1);
y = matrix(c(sort(1-exp(-x))),1);
i = matrix(c(1:N),1);
plot(y,i/N,type='s',col="red");
lines(matrix(c(0,1),1),matrix(c(0,1),1),col='blue');
dn_up = max(abs(y-i/N));
dn_down = max(abs(y-(i-1)/N));
dn = max(dn_up, dn_down);
j=matrix(-40:40,1);
p1 = 1-sum((-1)^j*exp(-2*(j*sqrt(N)*dn)^2));

p2=ks.test(y,punif)$p.value;

K = 1000;
DN=matrix(,1000)
for (k in 1:K)
{
	i = matrix(c(1:N),1);
	y=matrix(as.matrix(sort(matrix(c(matrix(c(runif(N)),1)),1))),1);
  DN[k] <- max(max(abs(y-i/N)), max(abs(y-(i-1)/N)));
}
p = sum(DN>=dn)/K;
p1
p2
p