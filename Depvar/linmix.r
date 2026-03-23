# linmix.r
source("gamrand.r");
##############<load a file>##################
file = as.matrix(read.csv("rats.csv"));     #
file = matrix(file[,],ncol=5);              #
d = c(151,199,246,283,320); #The first line #
rats = matrix(rbind(d,file),ncol=5);        #
#############################################
d = dim(rats)[1];
ni = dim(rats)[2];
n = d*ni;
y = matrix(t(rats),,1);
nloop = 11000;
burnin = 1000;
# storage
store_beta = matrix(rep(0,(nloop-burnin)*2),,2);
store_alpha = matrix(rep(0,(nloop-burnin)*d),,d);
store_var = matrix(rep(0,(nloop-burnin)*2),,2);
# priors
beta0 = matrix(c(100,10));
invVbeta = c(1/100, 1/100);
gamma0 = rbind(beta0,matrix(rep(0,d),,1));
nu_alpha = 3; lam_alpha = 100;
nu = 3; lam = 100;
# initialize the Markov chain
sigma2 = 100;
sigma2_alpha = 100;
# compute a few things before the loop
Z=kronecker(diag(d),matrix(rep(1,ni),,1));
xi = matrix(c(8, 15, 22, 29, 36));
X = kronecker(matrix(rep(1,d),,1), cbind(matrix(rep(1,ni),,1),xi));
W = cbind(X,Z);
WW = t(W) %*% W;
Wy = t(W) %*% y;
newnu_alpha = d/2 + nu_alpha;
newnu = n/2 + nu;
for (loop in 1:nloop){
  # sample alpha and beta
  invVgamma = diag(d+2)*c(invVbeta, 1/sigma2_alpha*(rep(1,d)));
  invDgamma = invVgamma + WW/sigma2;
  gammahat = solve(invDgamma, (invVgamma %*% gamma0 + Wy/sigma2));
  gamma = gammahat + solve(chol(invDgamma), matrix(rnorm(d+2)));
  beta = matrix(gamma[1:2]);
  alpha = matrix(gamma[3:nrow(gamma)]);
        # sample sigma2_alpha
  newlam_alpha = lam_alpha + sum(alpha^2)/2;
  sigma2_alpha = 1/gamrand(newnu_alpha, newlam_alpha);
        # sample sigma2
  newlam = lam + sum((y-W %*% gamma)^2)/2;
  sigma2 = 1/gamrand(newnu, newlam);
        # storage
  if (loop>burnin){
    i = loop-burnin;
    store_beta[i,] = t(beta);
    store_alpha[i,] = t(alpha);
    store_var[i,] = matrix(c(sigma2, sigma2_alpha),1);    
  }
}
betahat = matrix(rowMeans(t(store_beta)));
alphahat = matrix(rowMeans(t(store_alpha)));
varhat = matrix(rowMeans(t(store_var)));
print(betahat);
print(alphahat);
print(varhat);