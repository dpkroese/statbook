using Plots, LinearAlgebra, StatsBase,Distributions
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,3.8,4.8,1.7,-0.4,4.5,1.3,0.4,2.6,4,2.9,1.6]
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,8.05,3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04]
n = length(x)
b = 1.5

k(x,u,b) =  10*exp(-0.5*norm(x- u)^2/b^2) # scaled kernel
K = zeros(n,n)
for i = 1:n
   for j = 1:n
       K[i, j] = k(x[i], x[j], b)
   end
end

sigma=1 #starting value
sig2 = sigma^2
R = 10000 #number of samples for Gibbs sampler
sig2_store = zeros(R);
g_store = zeros(R,n);

for i=1:R  #gibbs sampler
   global sig2
    H = inv(K+sig2*I)
    mug = K'*H*y
    Sigg = K - K'*H*K
    U,S,V = svd(Sigg)
    gvec = mug + U*diagm(sqrt.(S))*randn(n)
    g_store[i,:] = gvec'
    sig2 = 1/rand(Gamma(n/2, 2/norm(y - gvec)^2))
    sig2_store[i] = sig2
end
println(sqrt(mean(sig2_store)))
# histogram(sqrt.(sig2_store))

