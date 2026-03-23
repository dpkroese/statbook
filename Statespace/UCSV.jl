using SparseArrays, Kronecker,LinearAlgebra, Distributions
using StatsBase,Plots,DelimitedFiles

function normpdf(Y,MU,SIG)
n, p = size(Y)
X = zeros(n,p)
for i=1:n
  for j=i:p
    X[i,j] = pdf(Normal(MU[i,j],SIG[i,j]),Y[i,j])
  end
end
return X
end


function SVRW(ystar,h,omega2h,Vh)
   T = length(h)
    # parameters for the Gaussian mixture
   pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575]
   mui = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]
   sig2i =[5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261]
   sigi = sqrt.(sig2i)
   s = zeros(Int64,T)
   for t=1:T
     q = zeros(7)
      for i=1:7
       q[i] = pi[i]*pdf(Normal(mui[i]-1.2704+h[t],sigi[i]),ystar[t])
      end
       q = q ./ sum(q)
       s[t] =  minimum(findall(cumsum(q) .> rand()))
   end
   
   H = sparse(I,T,T) - sparse(2:T,1:(T-1),ones(T-1),T,T)
   invOmegah = spdiagm(vec([1/Vh; 1/omega2h*ones(T-1,1)]))
   d = mui[s] 
   invSigystar = spdiagm(vec(1 ./ sig2i[s]))
   Kh = H'*invOmegah*H + invSigystar;
   R = cholesky(Kh)
   P = sparse(1:T,R.p,ones(T))
   Ch = P'*sparse(R.L)  
   hhat = Kh\(invSigystar*(ystar-d))
   h = hhat + Ch'\randn(T,1);
   return h, s
end

y = readdlm("USCPI.csv")
T = length(y)
nloop = 1000
Vtau = 9
Vh = 9
atau = 10
ltau = .25^2*(atau-1)
ah = 10 
lh = .2^2*(ah-1)
  # initialize the Markov chain
omega2tau = .25^2
omega2h = .2^2
h = log(var(y)*.8)*ones(T)
H = sparse(I,T,T) - sparse(2:T,1:(T-1),ones(T-1),T,T)
  # initialize for storage
store_omega2tau = zeros(nloop)
store_omega2h = zeros(nloop)
store_tau = zeros(nloop,T)
store_h = zeros(nloop,T)
 #  compute a few things
newatau = (T-1)/2 + atau
newah = (T-1)/2 + ah

for loop = 1:nloop
  global h,omega2tau,omega2h
invOmegatau = sparse(1:T,1:T, vec([1/Vtau 1/omega2tau*ones(1,T-1)]))
invSigy = sparse(1:T,1:T,vec(exp.(-h)))
Ktau = H'*invOmegatau*H + invSigy
R = cholesky(Ktau) # sparse Cholesky )
P = sparse(1:T,R.p,ones(T))
Ctau = P'*sparse(R.L) # Ctau*Ctau' = K
tauhat = Ktau\(invSigy*y)
tau = tauhat + Ctau'\randn(T)
ystar = log.((y-tau).^2 .+ .0001 )
h, s = SVRW(ystar,h,omega2h,Vh)
 # sample omega2tau
newltau = ltau + sum((tau[2:end]-tau[1:end-1]).^2)/2
omega2tau = 1/rand(Gamma(newatau, 1/newltau));
    # sample omega2h
newlh = lh + sum((h[2:end]-h[1:end-1]).^2)/2
omega2h = 1/rand(Gamma(newah, 1/newlh))
    store_tau[loop,:] = tau';
    store_h[loop,:] = h';
    store_omega2tau[loop] = omega2tau;
    store_omega2h[loop] = omega2h;
end

tauhat = mean(store_tau,dims=1)'
hhat = mean(store_h,dims=1)'

t = 1947.25:.25:2011.5
p1 = plot(1:length(tauhat),tauhat);
p2 = plot(1:length(hhat),hhat);
plot(p1,p2,)

# # testing
# T = length(h)
# # parameters for the Gaussian mixture
# pi = [0.0073 0.10556 0.00002 0.04395 0.34001 0.24566 0.2575]
# mui = [-10.12999 -3.97281 -8.56686 2.77786 0.61942 1.79518 -1.08819] .- 1.2704
# sig2i = [5.79596 2.61369 5.17950 0.16735 0.64009 0.34023 1.26261]
# sigi = sqrt.(sig2i)
# invOmegatau = sparse(1:T, 1:T, vec([1 / Vtau 1 / omega2tau * ones(1, T - 1)]))
# invSigy = sparse(1:T, 1:T, vec(exp.(-h)))
# Ktau = H' * invOmegatau * H + invSigy
# R = cholesky(Ktau) # sparse Cholesky )
# P = sparse(1:T, R.p, ones(T))
# Ctau = P' * sparse(R.L) # Ctau*Ctau' = K
# tauhat = Ktau \ (invSigy * y)
# tau = tauhat + Ctau' \ randn(T)
# ystar = log.((y - tau) .^ 2 .+ 0.0001)

# # sample s from a 7-point distrete distribution
# temprand = rand(T)

# q = zeros(7)
# for t=1:T
#   global q
#    for i=1:7
#     q[i] = pi[i]*pdf(Normal(mui[i]-1.2704,sigi[i]),ystar[t])
#    end
#     q = q ./ sum(q)
#     s[t] =  minimum(findall(cumsum(q) .> rand()))
# end
# H = sparse(I,T,T) - sparse(2:T,1:(T-1),ones(T-1),T,T)
#    invOmegah = spdiagm(vec([1/Vh; 1/omega2h*ones(T-1,1)]))
#    d = mui[s] 
#    invSigystar = spdiagm(vec(1 ./ sig2i[s]))
#    Kh = H'*invOmegah*H + invSigystar;
#    R = cholesky(Kh)
#    P = sparse(1:T,R.p,ones(T))
#    Ch = P'*sparse(R.L)  
#    hhat = Kh\(invSigystar*(ystar-d))
#    h = hhat + Ch'\randn(T,1);
#    return h, s