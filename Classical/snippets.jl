#Exercise
using StatsBase
N = 1000;
#comment one of the following pairs:
X = 4*rand(N) .- 2;  #method 1
H = 4*exp.(-X.^2/2); # method 1
#X = randn(N); # method 2
#H = sqrt(2*pi)*(abs.(X) .< 2); #  method 2
ell = mean(H)
sig = std(H)
RE = sig/sqrt(N)/ell
CI = ell*[1 - 1.96*RE 1 + 1.96*RE]

#Exercise
using Plots
n=10;
T1 = x ->  x/n - 1.96*sqrt((x/n)*(1-(x/n))/n);
T2 = x->  x/n + 1.96*sqrt((x/n)*(1-(x/n))/n);
p = 0:0.01:1;
cvp = zeros(size(p));
for i=1:length(p);
  tot = 0;
  for x = 0:n
      tot = tot + (T1(x) <= p[i] <= T2(x))*
      binomial(n,x)*p[i]^x*(1-p[i])^(n-x);
  end
  cvp[i] = tot;
end
plot(p,cvp)

#Exercise
using Distributions, Random, Statistics, LinearAlgebra, Plots
N = 100; x =  collect(1:N)/N ;
beta = [6; 13]; sigma = 2;# parameters
X = [ones(N,1) x];# design matrix
y = X*beta + sigma*randn(N,1); #draw the y-data
scatter(x,y); # plot the data
betahat = X'*X\(X'*y)  # solve the normal equations
sigmahat = norm(y - X*betahat)/sqrt(N) #estimate for sigma


# Q5.23
using LinearAlgebra, Statistics, Distributions
yy = [9 11 11.5 13.25;
9.5 10 12 11.5; 
9.75 10 9 12;
10 11.75 11.5 13.5;
13 10.5 13.25 11.5;
9.5 15 13 11.5]

n = length(yy); (nrow,ncol) = size(yy); y = vec(yy)
X_1 = ones(n,1)
KM = kron(diagm(ones(ncol)),ones(nrow,1)); X_2 = KM[:,1:ncol-1]
X_2[n-nrow+1:n,:] = -ones(nrow,ncol-1)
X = [X_1 X_2]
m = size(X,2);
betahat = X'*X\(X'*y)
ym = X*betahat
yk = X_1*mean(y); # omitting treatment effect
k=1
#number of parameters in reduced model
T = (n-m)/(m-k)*(norm(ym - yk)^2)/norm(y-ym)^2
pval = 1 - cdf(FDist(m-k,n-m),T)



# Q5.23
tquant = quantile(TDist(n-m),0.975) # 0.975 quantile
a = [1, 1, 0, 0]
ucl = a'*betahat + tquant*norm(y - X*betahat)*sqrt(a'*inv(X'*X)*a)/sqrt(n-m)
lcl = a'*betahat .- tquant*norm(y - X*betahat)*sqrt(a'*inv(X'*X)*a)/sqrt(n-m)
[lcl ucl]

# Q5.24
using LinearAlgebra, Statistics, Distributions
yy =  [13.8 11.7  14.0  12.6;
 12.9 16.7 15.5 13.8;
 25.9 29.8 27.8 25.0;
 18.0 23.1 23.0 16.9;
 15.2 20.2 19.9 13.7]

n = length(yy); (nrow,ncol) = size(yy); y = vec(yy)
X_1 = ones(n,1)
KM = kron(diagm(ones(ncol)),ones(nrow,1)); X_2 = KM[:,1:ncol-1]
X_2[n-nrow+1:n,:] = -ones(nrow,ncol-1)
X = [X_1 X_2]
m = size(X,2);
betahat = X'*X\(X'*y)
ym = X*betahat
yk = X_1*mean(y); # omitting treatment effect
k=1
#number of parameters in reduced model
T = (n-m)/(m-k)*(norm(ym - yk)^2)/norm(y-ym)^2
pval = 1 - cdf(FDist(m-k,n-m),T)
C = vcat(diagm(ones(nrow-1)), -ones(1,nrow-1))
X_3 = repeat(C,ncol,1)
X = [X_1 X_2 X_3]
m = size(X)[2]; #number of parameters in full model
betahat = X'*X\(X'*y) #estimate under the full model
ym = X*betahat
X_12 = [X_1 X_2] # omitting the block effect
k = size(X_12)[2] # number of parameters in reduced model
betahat_12 = X_12'*X_12\(X_12'*y)
y_12 = X_12*betahat_12;
T_12=(n-m)/(m-k)*(norm(y-y_12)^2 - norm(y-ym)^2)/norm(y-ym)^2
pval_12 = 1 - cdf(FDist(m-k,n-m),T_12)
X_13 = [X_1 X_3]; #omitting the treatment effect
k = size(X_13)[2]; #number of parameters in reduced model
betahat_13 = X_13'*X_13\(X_13'*y)
y_13 = X_13*betahat_13
T_13=(n-m)/(m-k)*(norm(y-y_13)^2 - norm(y-ym)^2)/norm(y-ym)^2
pval_13 = 1 - cdf(FDist(m-k,n-m),T_13)

