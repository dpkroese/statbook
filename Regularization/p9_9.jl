using LinearAlgebra,Statistics, Plots,LaTeXStrings
beta = [0.1, 0.2, 0.3]
X = [8.46  4.73   6.26;
     5.29  0.98  5.54;
     6.99  2.96  0.38;
     9.87  4.67  8.89;
     9.58  9.22  5.98]
r = [-0.83, 0.93, -0.24, -0.40, -0.02]
y = X*beta + r;

Xmu=mean(X,dims=1); ymu=mean(y);
y=y .- ymu; X=X .- Xmu;  # centering 

function CEmin(f, mu, sigma, N, Nel, tol)  # minimize function f via the CE method
n = length(mu) # dimension
   while maximum(sigma) > tol 
       ds = n == 1 ? sigma : diagm(sigma)  # scalars and vectors are treated differently
       X = randn(N,n)*ds .+ mu'            # generate N samples from N(mu,diag(sigma^2))
       fX= n==1 ? f.(X) : f.(eachrow(X));         # Compute the function values
       sortfX = sortslices(hcat(X, fX), dims=1, by = x -> x[n+1]) # sort rows by function val
       Elite = sortfX[1:Nel, 1:n]; # smallest (= elite) samples
       mu = n == 1 ? mean(Elite) : vec(mean(Elite,dims=1)) # take mean row-wise 
       sigma = n == 1 ? std(Elite) : vec(std(Elite,dims=1))     # take std row-wise
    end
    return f(mu), mu, sigma
end

function S(y,X,beta,lam)
    norm(y - X*beta)^2 + lam*norm(beta,1)
end

mu = [0,0,0]; sigma = 4.0*ones(3); N = 1000; Nel = 100; tol = 1E-5;

_,mu,_ = CEmin(beta-> S(y,X,beta,10),mu,sigma,N,Nel,tol)
print(mu)

n,p = size(X)
 lams = 0:0.1:30
 bs = zeros(p, length(lams))
 norms = zeros(length(lams))
 for i = 1:length(lams)
   mu = [0,0,0]; sigma = 4.0*ones(3); N = 1000; Nel = 100; tol = 1E-5;
   local minS,tol,Nel,N,sigma,mu
    lambda = lams[i]
    minS, mu, sigma = CEmin(beta -> S(y,X,beta,lambda),mu,sigma,N,Nel,tol)
    b = mu
    bs[:, i] = b
   norms[i] = norm(b, 1)
 end
 plot(tickfontsize=15,guidefontsize=20,tickfont = "Computer Modern")
 xlabel!(L"\mathrm{1-norm}")
 ylabel!(L"\bf{\beta}")
 for pp=1:p
      display(plot!(norms,bs[pp,:],legend=false,color="blue",linewidth=2))
 end
  plot!([0,1.3],[0,0],color=:black)