# probit_mle.r
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
Xbetat = X %*% betat;
phi = dnorm(Xbetat);  Phi = pnorm(Xbetat);
S = matrix(colSums(matrix(rep((y*phi/Phi-(1-y)*phi/(1-Phi)),each=k),byrow=TRUE,ncol=k)*X),k);
d = phi^2/(Phi*(1-Phi));
d=diag(1,n)*matrix(rep(d,each=n),byrow=TRUE,ncol=n);#make a sparse matrix
I=t(X)%*%(d)%*%X;#information matrix
betat = betat + solve(I,S);    
}
betat