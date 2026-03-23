# loglike_MA1.r
loglike_MA1<-function(theta,y){
  psi = theta[1]; sigma2 = theta[2]; 
  T = length(y);
  ####<Making a sparse matrix>####
  Sp=diag(T);                    #
  d=matrix(rep(0,T),1);          #
  sparse=rbind(d,Sp);            #
  sparse<-sparse[-(T+1),]        #
  ################################
  H = diag(T) + psi*sparse;
  HH = H  %*% t(H);
  l = (-T/2)*log(2*pi*sigma2) - (.5/sigma2) * t(y) %*% solve(HH,y);
  return (l);
}