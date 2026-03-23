rm(list=ls())
mu=3; sig=0.5;
alpha=0.05;n=10;
mu_count=0;sig_count=0;
for (k in 1:100)
{
  x=mu+matrix(rnorm(10),,1)*sig;
  mu_est = mean(x);
  sig_est = sd(as.numeric(x));
  tq = qt(1-alpha/2,n-1);
  mu_lo = mu_est - tq*sig_est/sqrt(n);
  mu_hi = mu_est +  tq*sig_est/sqrt(n);
  cq1=qchisq(1-alpha/2,n-1);
  cq2=qchisq(alpha/2,n-1);
  sig_lo = (n-1)*sig_est^2/cq1;
  sig_hi = (n-1)*sig_est^2/cq2;
  mu_count = mu_count + (mu > mu_lo & mu < mu_hi);
  sig_count =  sig_count + (sig^2 > sig_lo & sig^2 < sig_hi);
}
print(matrix(c(mu_count,sig_count),1,))