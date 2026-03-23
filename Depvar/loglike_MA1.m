% loglike_MA1.m
function l = loglike_MA1(theta,y)
psi = theta(1); sigma2 = theta(2); 
T = length(y);
H = speye(T) + psi*sparse(2:T,1:T-1,ones(1,T-1),T,T);
HH = H*H';
l = -T/2*log(2*pi*sigma2) - .5/sigma2*y'*(HH\y);
