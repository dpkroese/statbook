using DelimitedFiles, Distributions, LinearAlgebra
affair = readdlm("affair.csv", ',')
y = affair[:,1];
X = affair[:,2:end];
n, k = size(X);
#find the MLE and the information matrix
S = ones(k,1); #score
betat = (X'*X)\(X'*y); #initial guess
e = 10^(-5);  #tolerance level
while sum(abs.(S)) > e  #stopping criterion
    Xbetat = X*betat;
    phi = pdf.(Normal(0,1),Xbetat);
    Phi = cdf.(Normal(0,1),Xbetat);
    global S = sum(repeat(y.*phi./Phi-(1 .- y).*phi./(1 .-Phi),outer = [1,k]).*X,dims=1)';
    d = phi.^2 ./ (Phi.*(1 .-Phi));
    global IM = X'*diagm(vec(d))*X; # information matrix
    global betat = betat + IM\S;
end

burnin = 500;
nloop = 5000+burnin;
V = IM\I;  #scale matrix for the proposal
B = cholesky(Hermitian(V)).L;
nu = 5; #df for the proposal
b0 = zeros(k,1); #prior mean
V0 = 10*I; #prior covariance
# log-posterior density
logf(b)= (y'*log.(cdf.(Normal(0,1),X*b)) + (1 .-y)'*log.(cdf.(Normal(0,1),-X*b))- 0.5*(b-b0)'*(V0\(b-b0)))[1];
# log-proposal density
logprop(b)  = (-0.5*(k+nu)*log.(1 .+ (b-betat)'*(V\(b-betat))/nu))[1];
store_beta = zeros(nloop,k);
beta = betat;
for i = 1:nloop
# candidate draw from the t proposal
    global beta
    betac = betat + B*randn(k,1)*sqrt(nu/rand(Gamma(nu/2,2)));
    rho = logf(betac)-logf(beta) + logprop(beta)-logprop(betac);
    if exp(rho) > rand()
        beta = betac;
    end
    store_beta[i,:] = beta';
end
store_beta = store_beta[burnin+1:end,:]; # discard the burn-in
println(mean(store_beta,dims=1))
println(std(store_beta,dims=1))

# Marginal effects
N = size(store_beta,1);
store_ME = zeros(nloop-burnin,6);
xbar = mean(X,dims=1)';
for loop in 1:N
    global beta = store_beta[loop,:]; #ME for continuous variables
    store_ME[loop,[2 5]] = pdf.(Normal(0,1),xbar'*beta) .*beta[[3,6]];
    for j in [1 3 4 6] # ME for discrete variables 
        z0 = copy(xbar); z0[j+1] = 0;  # must use copy!
        z1 = copy(xbar); z1[j+1] = 1;
        store_ME[loop,j] = (cdf(Normal(0,1),z1'*beta) - cdf(Normal(0,1),z0'*beta))[1];
    end
end
mean(store_ME,dims=1)