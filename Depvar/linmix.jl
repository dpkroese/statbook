using DelimitedFiles, Kronecker, LinearAlgebra,SparseArrays, Distributions
rats = readdlm("rats.csv",',')
d, ni = size(rats)
n = d*ni
y = reshape(rats',n,1);
nloop = 11000
burnin = 1000
    # storage
store_beta = zeros(nloop-burnin,2)
store_alpha = zeros(nloop-burnin,d)
store_var = zeros(nloop-burnin,2)
    # priors
beta0 = [100 10]'
invVbeta = [1/100 1/100]
gamma0 = [beta0; spzeros(d,1)]
nu_alpha = 3; lam_alpha = 100
nu = 3; lam = 100;
    # initialize the Markov chain
sigma2 = 100
sigma2_alpha = 100
    # compute a few things before the loop
Z = sparse(kronecker(diagm(ones(d)), ones(ni,1)))
xi = [8 15 22 29 36]'
X = kronecker(ones(d,1),[ones(ni,1) xi])
W = [X Z]
WW = W'*W
Wy = W'*y
newnu_alpha = d/2 + nu_alpha;
newnu = n/2 + nu;
for loop=1:nloop
   global sigma2_alpha,sigma2
        # sample alpha and beta
    invVgamma = sparse(1:d+2,1:d+2, vec([invVbeta 1/sigma2_alpha*ones(1,d)]))
    invDgamma = invVgamma + WW/sigma2
    gammahat = invDgamma\(invVgamma*gamma0 + Wy/sigma2)
    gamma = gammahat + cholesky(invDgamma).L'\randn(d+2,1)
    beta = gamma[1:2]
    alpha = gamma[3:end]
        # sample sigma2_alpha
    newlam_alpha = lam_alpha + sum(alpha.^2)/2
    sigma2_alpha = 1/rand(Gamma(newnu_alpha, 1/newlam_alpha));
        # sample sigma2
    newlam = lam + sum((y-W*gamma).^2)/2
    sigma2 = 1/rand(Gamma(newnu, 1/newlam))
        # storage
    if loop>burnin
        i = loop-burnin;
        store_beta[i,:] = beta';
        store_alpha[i,:] = alpha';
        store_var[i,:] = [sigma2 sigma2_alpha];    
    end
end
betahat = mean(store_beta,dims=1)
alphahat = mean(store_alpha,dims=1)
varhat = mean(store_var,dims=1)
