% lmvnpdf.m
function lden = lmvnpdf(x,mu,Sig)
m = size(mu,1);
e = x-mu;
c = -m/2*log(2*pi) - sum(log(diag(chol(Sig,'lower'))));
lden = c - .5*e'*(Sig\e);