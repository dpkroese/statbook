# probit_bayes_gibbs.r
source("tnormrnd.r")
affair = as.matrix(read.csv("affair.csv"));
d = c(0,1,1,10,0,0,18,1); #The first line
affair = rbind(d,affair);
#if you get an error like 'No such file or directory' 
#need to set working directory before run this.
y = matrix(affair[,1],,1);
X = matrix(affair[,(2:ncol(affair))],,7);
XX = t(X) %*% X;
n = dim(X)[1];
k = dim(X)[2];
V0 = 10*diag(1,k);            # prior covariance
invV0 = solve(V0,diag(1,k));
burnin = 500;
nloop = 5000+burnin;
store_beta = matrix(rep(0,nloop*k),ncol=k)
z = y;                     # initial guess
beta = solve(XX,(t(X) %*% z));
                           # compute a few things before the loop
id0 = which(y==0); 
id1 = which(y==1);
n0 = length(id0);
n1=n-n0;
V = solve((invV0 + XX),diag(1,k)); # posterior covariance
for (i in 1:nloop)
{           # sample y^*
  Xb = X %*% beta;
  z[id0] = tnormrnd(Xb[id0],matrix(rep(1,n0)),-Inf,0);  
  z[id1] = tnormrnd(Xb[id1],matrix(rep(1,n1)),0,Inf); 
            # sample beta
  dbeta = t(X) %*% z;
  beta = V %*% dbeta + t(chol(V)) %*% matrix(rnorm(k));
  store_beta[i,] = t(beta);
}
store_beta = store_beta[((burnin+1):nrow(store_beta)),];  # discard the burnin
colMeans(store_beta)
apply(store_beta,2,sd)