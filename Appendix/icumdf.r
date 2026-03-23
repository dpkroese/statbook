icumdf<- function(dist, y, varargin){
  switch(dist,
         'norm'={
           mu <- varargin[[1]];
           sigma <- varargin[[2]];
           x <- (mu  + sigma*sqrt(2)*qnorm((1+(2*y-1))/2)/sqrt(2));   
         },
         't'={
           nu <- varargin[[1]];
           x <- sqrt(nu/qbeta(2*(1-y),nu/2,1/2) - nu);
         },
         'gamma'={
           alpha <- varargin[[1]];
           lambda <- varargin[[2]];
           # different from Stats toolbox
           x <- qgamma(y,alpha)/lambda; 
          },
         'chi2'={
           n <- varargin[[1]];
           x <- qgamma(y,n/2)*2;
         },
         'F'={
           m <- varargin[[1]];
           n <- varargin[[2]];
           x <- n/m/qbeta(1-y,n/2,m/2) - n/m;
         }
  )
  return (x);
}