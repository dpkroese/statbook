# bioassay_mle.r
y = matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,1,1,1,1));
x = matrix(rep(c(-0.863,-0.296,-0.053,0.727),each=5),ncol=4);
X = cbind(matrix(rep(1,20),,1),c(x));      # design matrix
betat = solve((t(X)%*%X),(t(X)%*%y));      # initial guess
S = matrix(rep(1,2),,1);    # score
e = 10^(-5);                # tolerance level
while (sum(abs(S)) > 10^(-5))# stopping criteria
{
  mu = 1/(1+exp(-X%*%betat));
  S = matrix(colSums(cbind(y-mu,y-mu)*X),,1);
  d=diag(1,20)*matrix(rep(mu*(1-mu),each=20),byrow=TRUE,ncol=20);#make a sparse matrix
  I=t(X)%*%(d)%*%X;#information matrix
  betat = betat + solve(I,S);   
} 
V = solve(I,diag(1,2));               # covariance matrix            
print(V)