function loglike_ARMA11(psi,X,y)
  T = length(y)
  H = sparse(I,T,T) + psi*sparse(2:T,1:T-1,ones(T-1),T,T)
  HH = H*H'
  rhohat = (X'*(HH\X))\(X'*(HH\y))
  uhat = y-X*rhohat
  sigma2hat = uhat'*(HH\uhat)/T
  l = -T/2*log(2*pi*sigma2hat) - .5/sigma2hat*uhat'*(HH\uhat)
  return l, rhohat, sigma2hat
end
