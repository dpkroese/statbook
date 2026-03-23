using SparseArrays, Kronecker,LinearAlgebra, Distributions
using StatsBase,Plots,DelimitedFiles

function SURform(X)  # same as matlab
   r, c = size(X);
   idi = vec(Int64.(kronecker(1:r,ones(c))))
   idj = 1:r*c
   return sparse(idi,idj,X'[:])  #maybe X?
end

function gamrnd(a,c) #c is scale not rate
  n = length(c)
  x = zeros(n)
  for i=1:n
   x[i]  = rand(Gamma(a,c[i]))
  end
  return x
end

USCPI= readdlm("USCPI.csv")
nloop = 10000
p = 2;    # number of lags
y0 = USCPI[1:p]
y = USCPI[p+1:end] 
T = length(y)
q = p+1  # dimension of states
Tq = T*q 
    # priors
asigma2 = 5 
lsigma2 = 1*(asigma2-1)
aomega2 = 5
lomega2 = (aomega2-1)*[0.5^2; 0.1^2*ones(p,1)]  
invOmega0 = ones(q)/5
   # initialize
omega2 = .1*ones(q)
sigma2 = 1
store_omega2 = zeros(nloop,q)
store_sigma2 = zeros(nloop)
store_beta = zeros(Tq)
    # construct/compute a few things
X = [ones(T,1) [y0[end]; y[1:end-1]] [y0; y[1:end-2]]]
bigX = SURform(X)
H = sparse(I,Tq,Tq) - sparse(q+1:Tq,1:(T-1)*q, ones((T-1)*q),Tq,Tq)
newaomega2 = aomega2 + T - 1
newasigma2 = asigma2 + T 
for loop = 1:nloop
   global omega2,sigma2,store_beta,betahat
        # sample beta
    invS = sparse(1:Tq,1:Tq,vec([invOmega0' repeat(1 ./ omega2',1,T-1)]))
    K = H'*invS*H + bigX'*bigX/sigma2
    R = cholesky(K)  # sparse Cholesky 
    P = sparse(1:Tq,R.p,ones(Tq))
    C = P'*sparse(R.L) # C*C' = K
    betahat = K\(bigX'*y/sigma2)
    beta = betahat + C'\ randn(Tq)    
       # sample omega2
    erromega2 = reshape(H*beta,q,T)
    newlomega2 = lomega2 + sum(erromega2[:,2:end].^2,dims=2)/2
    omega2 = 1 ./ gamrnd(newaomega2, 1 ./ newlomega2)    
        # sample sigma2
    newlsigma2 = lsigma2 + sum((y-bigX*beta).^2)/2
    sigma2 = 1/rand(Gamma(newasigma2,1/newlsigma2))   
        # store 
    store_beta = store_beta + beta 
    store_omega2[loop,:] = omega2'
    store_sigma2[loop] = sigma2 
end
betahat = store_beta/nloop
sigma2hat = mean(store_sigma2)
omega2hat = mean(store_omega2,dims=1)
t = 1947.25:.25:2011.5
p1 = plot(t[3:end],betahat[1:3:end]);
p2 = plot(t[3:end],betahat[2:3:end]);
p3 = plot(t[3:end],betahat[3:3:end]);
plot(p1,p2,p3,layout=(1,3))