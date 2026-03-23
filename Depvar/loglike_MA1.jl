function loglike_MA1(theta,y)
  # the log-likelihood function for MA(1)
  # input: theta = [psi sigma2]; y = data
  psi = theta[1]; sigma2 = theta[2];
  T = length(y);
  H = sparse(I,T,T) .+ psi*sparse(2:T,1:T-1,ones(T-1),T,T);
  HH = H*H';
  l = -T/2*log(2*pi*sigma2) - .5/sigma2*y'*(HH\y);
  return l
end
