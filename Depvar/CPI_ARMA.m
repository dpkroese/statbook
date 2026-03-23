% CPI_ARMA.m
load USCPI.csv
y0 = USCPI(1);
y = USCPI(2:end);
T = length(y);
X = [ones(T,1) [y0; y(1:end-1)]];
f = @(psi) -loglike_ARMA11(psi,X,y);
psihat = fminsearch(f,0);
[l rhohat sig2hat] = loglike_ARMA11(psihat,X,y)

