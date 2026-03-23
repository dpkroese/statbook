function AR1_loglike(rho,y,X)
T = length(y);
H = sparse(I,T,T) .- rho*sparse(2:T,1:T-1,ones(T-1),T,T);
HH = H'*H;
betahat = (X'*HH*X)\(X'*HH*y);
e = y-X*betahat;
sigma2hat = e'*HH*e/T;
l = -T/2*log(2*pi*sigma2hat) - .5/sigma2hat*e'*HH*e;
return l, betahat,sigma2hat
end
