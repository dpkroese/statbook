loglike_ARMA11<-function(psi,X,y){
  T = length(y);
  ####<Making a sparse matrix>####
  Sp=diag(T);                    #
  d=matrix(rep(0,T),1);          #
  sparse=rbind(d,Sp);            #
  sparse<-sparse[-(T+1),]        #
  ################################
  H = diag(T) + psi*sparse;
  HH = H  %*% t(H);
  rhohat = solve((t(X) %*% solve(HH,X)),(t(X) %*% solve(HH,y)));
  uhat = y-X %*% rhohat;
  sig2hat = (t(uhat) %*% solve(HH,uhat)/T)[1];    
  l = -T/2*log(2*pi*sig2hat) - .5/sig2hat*t(uhat) %*% solve(HH,uhat);
  result=list(l,rhohat,sig2hat);
  return (result);
}