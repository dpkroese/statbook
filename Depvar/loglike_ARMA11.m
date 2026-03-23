% loglike_ARMA11.m
function [l rhohat sig2hat] = loglike_ARMA11(psi,X,y)
T = length(y);
H = speye(T) + psi*sparse(2:T,1:T-1,ones(1,T-1),T,T);
HH = H*H';
rhohat = (X'*(HH\X))\(X'*(HH\y));
uhat = y-X*rhohat;
sig2hat = uhat'*(HH\uhat)/T;    
l = -T/2*log(2*pi*sig2hat) - .5/sig2hat*uhat'*(HH\uhat);

