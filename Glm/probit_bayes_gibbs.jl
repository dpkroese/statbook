using DelimitedFiles, Distributions, LinearAlgebra
affair = readdlm("affair.csv", ',')
y = Int64.(affair[:,1]); # convert the response to integers
X = affair[:,2:end];
XX = X'*X;
n,  k = size(X);
V0 = 10*diagm(ones(k)); #prior covariance
invV0 = V0\I;
burnin = 500;
nloop = 5000+burnin;
store_beta = zeros(nloop,k);
z = 1.0*y; # initial guess, new float copy
beta = XX\(X'*z);
# compute a few things before the loop
id0 = findall(y .== 0); id1 = findall(y .==1);
n0 = length(id0); n1=n-n0;
V = (invV0 + XX)\I; # posterior covariance
for i in  1:nloop 
    # sample z
    global beta
    global Xb = X*beta;
    for k in id0
        z[k] = rand(Truncated(Normal(Xb[k],1), -Inf,0))
    end
    for k in id1
        z[k] = rand(Truncated(Normal(Xb[k],1), 0, Inf))
    end
    # sample beta
    dbeta = X'*z;
    beta = V*dbeta + cholesky(Hermitian(V)).L*randn(k,1);
    store_beta[i,:] = beta';
end
store_beta = store_beta[burnin+1:end,:]; # discard the burn-in
mean(store_beta,dims=1)

# for Problem 9.4
N = size(store_beta)[1];
store_prob = zeros(N);
x = [1 1 10 1 0 16 1];
for loop=1:N
  store_prob[loop] = cdf.(Normal(0,1),x*store_beta[loop,:])[1]
end
include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using .ThetaKDE
kde(store_prob,plt=true)