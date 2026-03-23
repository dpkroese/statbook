# probit_bayes.r
source("gamrand.r")
affair=as.matrix(read.csv("affair.csv"));
d = c(0,1,1,10,0,0,18,1); #The first line
affair = rbind(d,affair);
#if you get an error like 'No such file or directory' 
#need to set working directory before run this.
y = matrix(affair[,1],,1);
X = matrix(affair[,(2:ncol(affair))],,7);
n = dim(X)[1];
k = dim(X)[2];
## obtain the proposal density for the MH step
S = matrix(rep(1,2),,1);                   # score
betat = solve((t(X)%*%X),(t(X)%*%y));      # initial guess
e = 10^(-5);                               # tolerance level
while (sum(abs(S)) > e){        # stopping criteria
  Xbetat = X%*%betat;
  phi = dnorm(Xbetat);  Phi = pnorm(Xbetat);
  S = matrix(colSums(matrix(rep((y*phi/Phi-(1-y)*phi/(1-Phi)),each=k),byrow=TRUE,ncol=k)*X),k);
  d = phi^2/(Phi*(1-Phi));
  d=diag(1,n)*matrix(rep(d,each=n),byrow=TRUE,ncol=n);#make a sparse matrix
  I=t(X)%*%(d)%*%X;#information matrix
  betat = betat + solve(I,S);    
}
burnin = 500;
nloop = 5000+burnin;
V = solve(I,diag(1,k)) # scale matrix for the proposal
B = t(chol(V));
nu = 5;                     # df for the proposal
b0 = matrix(rep(0,k),ncol=1);            # prior mean
V0 = 10*diag(1,k);             # prior covariance
# log-posterior density
logf<-function(b){t(y)%*%log(pnorm(X%*%b)) + t(1-y)%*%log(pnorm(-X%*%b)) - .5*t(b-b0)%*%(solve(V0,(b-b0)))};
# log-proposal density
logprop<-function(b){-0.5*(k+nu)*log(1+t(b-betat)%*%(solve(V,(b-betat))/nu))} ;
store_beta=matrix(rep(0,nloop*k),ncol=k)
beta = betat;
for (i in 1:nloop)
{
  # candidate draw from the t proposal
  betac = betat + B%*%matrix(rnorm(k),,1)*sqrt(nu/gamrand(nu/2,1/2)); 
  rho = logf(betac)-logf(beta) +logprop(beta)-logprop(betac);
  if (exp(rho) > runif(1)){
    beta = betac;
  }
  store_beta[i,] = t(beta);    
}
store_beta = store_beta[((burnin+1):nrow(store_beta)),];  # discard the burnin
colMeans(store_beta)