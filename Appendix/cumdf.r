cumdf<- function(dist, x, varargin){
	switch(dist,
	'norm'={
		mu <- varargin[[1]]; 
    sigma <- varargin[[2]];
    f=((x - mu)/sigma )/sqrt(2);
    y <- (((2*pnorm(f*sqrt(2))-1)+1)/2);
	},
	't'={
		nu <- varargin[[1]];
		y <- (1-0.5*pbeta(nu/(nu+x^2),nu/2,1/2));
	},
	'gamma'={
		alpha <- varargin[[1]];
		lambda <- varargin[[2]];
		y <- pgamma(lambda*x,alpha); 
	},
  'chi2'={
    n <- varargin[[1]];
    y <- pgamma(x/2,n/2);  
  },
    'F'={
    m <- varargin[[1]];
    n <- varargin[[2]];
    y <- (1 - pbeta(n/(n+m*x),n/2,m/2));       
  }
)
  return (y);
}