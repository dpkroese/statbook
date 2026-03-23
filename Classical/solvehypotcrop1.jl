using LinearAlgebra, Statistics, Distributions
yy = [9.2988 9.4978 9.7604 10.1025;
 8.2111 8.3387 8.5018 8.1942; 
 9.0688 9.1284 9.3484 9.5086;
 8.2552 7.8999 8.4859 8.9485]
n = length(yy); (nrow,ncol) = size(yy); y = vec(yy)
X_1 = ones(n,1)
KM = kron(diagm(ones(ncol)),ones(nrow,1)); 
X_2 = KM[:,1:ncol-1]
X_2[n-nrow+1:n,:] = -ones(nrow,ncol-1)
X = [X_1 X_2]
m = size(X,2)
betahat = X'*X\(X'*y)
ym = X*betahat;
yk = X_1*mean(y); # omitting treatment effect
k=1;  #number of parameters in reduced model
T = (n-m)/(m-k)*(norm(ym - yk)^2)/norm(y-ym)^2
pval = 1 - cdf(FDist(m-k,n-m),T)