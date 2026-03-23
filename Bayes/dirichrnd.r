source("gamrand.r");
dirichrnd<-function(alpha){
  n=length(alpha)-1;
  Y=rep(NaN,n+1);
  for (k in (1:(n+1))){
    Y[k]=gamrand(alpha[k],1);  
  }
  x=Y[1:n]/sum(Y);
}