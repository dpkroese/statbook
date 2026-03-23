rm(list=ls());
xbar_obs = -0.7; s_obs = 0.4;
t_obs = 2*xbar_obs/s_obs;
N = 10^5;
tot = 0;
for (i in 1:N)
{
  x = matrix(c(rnorm(4)),,1);
  s = sd(as.numeric(x)); xbar=mean(x);
  t = 2*xbar/s;
  tot = tot + as.numeric(t <= -3.5);
}
phat = tot/N
phat
pt(-3.5,3)