kde<-function(data,n,MIN,MAX,len,loop) {
#version with "loop" variable for plotting
  data<-matrix(data,,1);
  if (nargs()<4){
    n=2^14;
  }
  n=2^ceiling(log2(n));
  if (nargs()<6){
    minimum=min(data);
    maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/10;
    MAX=maximum+Range/10;
  }
  R=MAX-MIN;
  dx=R/(n-1);
  density=0;
  cdf=0;
  xmesh=MIN+matrix(seq(0,R,dx),1);
  N=length(unique(data))
  initial_data=matrix((hist(data,xmesh,plot=FALSE)$counts)/N);
  initial_data=rbind(initial_data,0) # The last term counts the out of bin values
  initial_data=initial_data/sum(initial_data);
  a=dct1d(initial_data);
  assign("I", matrix(seq(1,n-1))^2)
  a2=matrix((a[2:nrow(a)]/2)^2);
  tryCatch(
    (t_star=uniroot(fixed_point,c(0,0.1),x=N,y=I,z=a2)$root),
    error = (t_star=.28*N^(-2/5))
  )
  a_t<-a*exp(-matrix(seq(0,n-1))^2*pi^2*t_star/2);
  if (len>1|len==0){
    density=idct1d(a_t)/R
  }
  bandwidth=sqrt(t_star)*R;
  if (len==0){
    if(loop==1){
      ##############################
      win.graph() #---for Windows##
      #quartz() #---for Mac       ###
      ##############################
      plot(smooth.spline(xmesh,density),type='l',col='blue',ylab='',xlab='')
    }
    else{
      lines(xmesh,density,col='blue')
     }
  }
  if (len>=3){
    f=2*pi^2*sum(I*a2*exp(-I*pi^2*t_star));
    t_cdf<-(sqrt(pi)*f*N)^(-2/3);
    a_cdf=a*exp(-matrix(0:(n-1))^2*pi^2*t_cdf/2);
    cdf<-matrix(cumsum(idct1d(a_cdf))*(dx/R));
    bandwidth_cdf=sqrt(t_cdf)*R;
  }
  result=list(bandwidth,density,xmesh,cdf);
  return (result)
}
#################################################################
fixed_point<-function(t,x,y,z) { 
  #x=N,y=I,z=a2
  l=7;
  f=2*pi^(2*l)*sum(y^l*z*exp(-y*pi^2*t));
  for(s in seq(l-1,2,-1)){
    K0=prod(seq(1,2*s-1,2))/sqrt(2*pi);
    const=(1+(1/2)^(s+1/2))/3;
    time=(2*const*K0/x/f)^(2/(3+2*s));
    f=2*pi^(2*s)*sum(y^s*z*exp(-y*pi^2*time));
  }
  return (t-(2*x*sqrt(pi)*f)^(-2/5));
}

idct1d<-function(data_1) {
  nrows=nrow(data_1);
  ncols=ncol(data_1);
  weights=t(nrows*exp(1i*matrix(seq(0,nrows-1),1)*pi/(2*nrows)));
  data_1=Re(fft(weights*data_1)/length(data_1))
  out=matrix(rep(0,nrows));
  out[seq(1,nrows,2)]=matrix(data_1[1:(nrows/2)],1);
  out[seq(2,nrows,2)]=matrix(data_1[seq(nrows,nrows/2+1,-1)],1);
  return (out);
}

dct1d<-function(data_1) {
  nrows=nrow(data_1);
  ncols=ncol(data_1);
  weight = rbind(1,2*t(exp(-1i*matrix(seq(1,nrows-1),1)*pi/(2*nrows))));
  data_1 = rbind(matrix(data_1[seq(1,nrows,2)]),matrix(data_1[seq(nrows,2,-2)]));
  data_1 = Re(weight*fft(data_1))
  return (data_1);  
}