rm(list=ls());
yy=rbind(c(9.2988, 9.4978, 9.7604, 10.1025),c(8.2111, 8.3387, 8.5018,  8.1942), c(9.0688, 9.1284, 9.3484,  9.5086),c(8.2552, 7.8999, 8.4859,  8.9485));

n=as.numeric(length(yy)); nrow<-as.numeric(dim(yy)[1]); ncol<-as.numeric(dim(yy)[2]);
y=matrix(as.numeric(yy),,1);
X_1 = matrix(1,n,1);
KM = kronecker(diag(ncol),matrix(1,nrow,1));
X_2 = KM[1:n,0:(ncol-1)];
X_2[(n-nrow+1):n,0:3]=-matrix(1,nrow,ncol-1);
C = t(cbind(diag(nrow-1),-matrix(1,nrow-1,1)));
#make a function
repmat<-function(a,n,m) {kronecker(matrix(1,n,m),a)};
X_3=repmat(C,ncol,1);
X=cbind(X_1,X_2,X_3);
m = as.numeric(dim(X)[2]);
betahat=solve(t(X)%*%(X))%*%(t(X)%*%y);
ym=X%*%betahat;
X_12=cbind(X_1,X_2);
k=as.numeric(dim(X_12)[2]);
betahat_12 = solve(t(X_12)%*%(X_12))%*%(t(X_12)%*%y);
y_12=X_12%*%betahat_12;
T_12 = (n-m)/(m-k)*(norm(y-y_12,"F")^2 - norm(y-ym,"F")^2)/norm(y-ym,"F")^2
pval_12 = 1-pf(T_12,m-k,n-m);
cat("pval_12=\n"); print(pval_12)
X_13=cbind(X_1,X_3);
k = as.numeric(dim(X_13)[2]);
betahat_13 = solve(t(X_13)%*%(X_13))%*%(t(X_13)%*%y);
y_13 = X_13%*%betahat_13;
T_13 = (n-m)/(m-k)*(norm(y-y_13,"F")^2 - norm(y-ym,"F")^2)/norm(y-ym,"F")^2
pval_13 = 1-pf(T_13,m-k,n-m)
cat("pval_13=\n"); print(pval_13)