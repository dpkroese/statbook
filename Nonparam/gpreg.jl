using Plots, LinearAlgebra, StatsBase
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,3.8,4.8,1.7,-0.4,4.5,1.3,0.4,2.6,4,2.9,1.6]
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,8.05,3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04]
n = length(x)

# construct the optimal cubic regression polynomial
X = ones(n)
for k=1:3
   global X
   X = hcat(X, x.^k)   # make the design matrix    
end
beta = X'*X\(X'*y)  # optimal parameters

g(x, beta) =  beta[1] .+ beta[2]*x .+ beta[3]*x.^2 .+ beta[4]*x.^3 

b = 1.75 
sigma = 0.84

k(x,u,b) =  10*exp(-0.5*norm(x- u)^2/b^2) # scaled kernel

function gp_pred(xtest, xtrain, ytrain, sigma, b, k)
    T = length(xtrain); N = length(xtest)
    K = zeros(T, T)
    mu = zeros(N); sigma_squared = zeros(N)
    for z = 1:N
        kappa = zeros(T)
        for i = 1:T
            kappa[i] = k(xtest[z], xtrain[i], b)
            for j = 1:T
                K[i, j] = k(xtrain[i], xtrain[j], b)
            end
        end
        mu[z] = kappa' * (K + sigma^2 * I)^(-1) * ytrain
        sigma_squared[z] = k(xtest[z], xtest[z], b) - kappa' * (K + sigma^2 * I)^(-1) * kappa
    end
    return mu, sigma_squared
end

xtilde = minimum(x):0.01:maximum(x)
mu,sigsquare = gp_pred(xtilde,x,y,sigma,b,k)
lo = mu -1.96*sqrt.(sigsquare); hi = mu+1.96*sqrt.(sigsquare);
plot(xtilde, mu,lw=2,color=:red,legend=false)
plot!(xtilde, mu, ribbon = (lo .- hi) ./ 2, fillalpha = 0.1,color=:black)
plot!(xtilde, mu,lw=2,color=:red)
plot!(xtilde,g(xtilde,beta),lw=2,color=:blue,ls=:dash)
scatter!(x,y,color=:black)
