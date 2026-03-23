invgampdf<-function(z,a,l){
  l^a*z^(-a-1)*exp(-l/z)/gamma(a);
}