rm(list=ls());
cat("\f")
graphics.off()
source("kde.r");
source("gamrand.r");
source("invgampdf.r");
n=10;
#set.seed(1234);
X=matrix(rnorm(n),1); #generate the data
sample_mean=mean(X);
sample_var= var(X[,]);
sig2=var(X[1:10]);
mu=sample_mean; #initial state
N=10^5;# Gibbs sample size
gibbs_sample=matrix(0,N,2);
for (k in 1:N){
  mu=sample_mean + sqrt(sig2/n)*rnorm(1);  #draw mu
  V = sum((X-mu)^2)/2;
  sig2=1/gamrand(n/2,V);
  gibbs_sample[k,]=cbind(mu,sig2);    
}    
# t=rbind(seq(-3,3,0.01));
# y = sqrt(n)/sd(X)*pt(sqrt(n)*(t-sample_mean)/sd(X),n-1)
# plot(t,y,col='red',type='l');
# lines(c(sample_mean,sample_mean),c(0,max(y)));
# kde(gibbs_sample[,1])
t = rbind(seq(0,5,0.01));
y = invgampdf((t/(n-1)/sample_var),((n-1)/2),1/2)/(n-1)/sample_var;
#plot(t,y,col='red')
result<-kde(gibbs_sample[,2],len=0,loop=1); # plot the histogram of the variance
lines(c(sample_var,sample_var),c(0,max(y,na.rm=TRUE)),col='blue');
#mean(gibbs_sample[,2])
#sample_var
