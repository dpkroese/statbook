x=matrix(c(125,18,20,34));
n=sum(x);
theta=4*(x[1]/n-1/2);
err=1;
while (abs(err) > 10^(-5))
{
  z=2*x[1]/(2+theta);
  temp=(x[1]+x[4]-z)/(n-z);
  err = theta - temp;
  theta=temp;
  print(theta);
}