using LinearAlgebra, Distributions
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,
3.8,4.8,1.7,-0.4,4.5,1.3,0.4,2.6,4,2.9,1.6];
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,8.05,
3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04];
n = length(x);
X = hcat(ones(n), x, x.^2, x.^3, x.^4, x.^5);
XX = X'*X;
Xy = X'*y;
m = 6;
N = 10^5;   # Gibbs sample size
IM = diagm(ones(m))         
V0 = 100*IM     # prior for beta
invV0 = V0\IM;
alp0 = 2; lam0 = 1;  # prior for sig2
beta = XX\Xy;
sig2 = sum((y -X*beta).^2)/n
gibbs_sample = zeros(N,m+1);
lpostden_sample = zeros(N,4); 

for k in 1:N   
   global beta, sig2 
    D = (invV0 + XX/sig2)\IM;
    betahat = D*(Xy/sig2)
    beta = betahat + cholesky(Hermitian(D)).L*randn(m);
    sig2 = 1/rand(Gamma(alp0+n/2,1/(lam0+sum((y-X*beta).^2)/2)));
    gibbs_sample[k,:]=[beta' sig2];
    lp1 = logpdf(MvNormal(betahat[3:end],
                  Hermitian(D[3:end,3:end])),zeros(4))
    lp2 = logpdf(MvNormal(betahat[4:end],
                  Hermitian(D[4:end,4:end])),zeros(3))
    lp3 = logpdf(MvNormal(betahat[5:end],
                  Hermitian(D[5:end,5:end])),zeros(2))
    lp4 = logpdf(Normal(betahat[6],sqrt(D[6,6])),0)
    lpostden_sample[k,:] = [lp1 lp2 lp3 lp4];
end
lpostden = zeros(4,1);
for i in 1:4
    maxpden = maximum(lpostden_sample[:,i]);
    lpostden[i]=log.(mean(exp.(lpostden_sample[:,i].-maxpden))) + maxpden;
end
lpriden = zeros(4,1);
lpriden[1] = logpdf(MvNormal(zeros(4),
                      V0[3:end,3:end]),zeros(4));
lpriden[2] = logpdf(MvNormal(zeros(3),
                      V0[4:end,4:end]),zeros(3));
lpriden[3] = logpdf(MvNormal(zeros(2),
                      V0[5:end,5:end]),zeros(2));
lpriden[4] = logpdf(Normal(0,sqrt(V0[6,6])),0);
lBF = lpostden - lpriden;