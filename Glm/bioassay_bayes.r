# bioassay_bayes.r
source("gamrand.r")
y = matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,1,1,1,1));
x = matrix(rep(c(-0.863,-0.296,-0.053,0.727),each=5),ncol=4);
A = cbind(matrix(rep(1,20),,1),c(x));      # design matrix
betat = solve((t(A)%*%A),(t(A)%*%y));      # initial guess
S = matrix(rep(1,2),,1);    # score
e = 10^(-5);                # tolerance level
while (sum(abs(S)) > e){  # stopping criteria
  mu = 1/(1+exp(-A%*%betat));
  S = matrix(colSums(cbind(y-mu,y-mu)*A),,1);
  I = t(A)%*%diag(c(mu*(1-mu)))%*%A; # infomation matrix
  betat = betat + solve(I,S);   
} 
V = solve(I,diag(1,2))
B = t(chol(V));
burnin = 500;
nloop = 10000+burnin; 
store_beta=matrix(rep(0,2*nloop),ncol=2)
nu = 5;  		     # df for the proposal 
# log posterior density
logf<-function(b){sum((y-1)*(A%*%b)-log((1+exp(-A%*%b))))};
# log density of the t proposal
logprop<-function(b){-0.5*(nu+2)*log(1+t(b-betat)%*%(solve(V,(b-betat))/nu))} ;
beta = betat;               # initialize the chain
for (i in 1:nloop)
{
  # candidate draw from the t proposal
  betac = betat + B%*%matrix(rnorm(2),,1)*sqrt(nu/gamrand(nu/2,1/2)); 
  rho = logf(betac)-logf(beta) +logprop(beta)-logprop(betac);
  if (exp(rho) > runif(1)){
    beta = betac;
  }
  store_beta[i,] = t(beta);    
}
store_beta = store_beta[((burnin+1):nrow(store_beta)),];  # discard the burnin
colMeans(store_beta)
cov(store_beta)