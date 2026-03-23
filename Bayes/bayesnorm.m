%bayesnorm.m
clear all, clc, clf
randn('state',1234);
rand('twister',3421);
n=10; 
X=randn(1,n); %generate the data
sample_mean=mean(X);
sample_var= var(X);
sig2=var(X); mu=sample_mean; %initial state
N=10^5;% Gibbs sample size
gibbs_sample=zeros(N,2);
for k=1:N    
    mu=sample_mean + sqrt(sig2/n)*randn;  %draw mu
    V = sum((X -mu).^2)/2;
    sig2=1/gamrand(n/2,V);
    gibbs_sample(k,:)=[mu,sig2];  
end
hold on
%t = [-3:0.01:3]
%y = sqrt(n)/std(X)*pdf('t',sqrt(n)*(t-sample_mean)/std(X),n-1)
%plot(t,y,'r')
%line([sample_mean, sample_mean],[0,max(y)])
%kde(gibbs_sample(:,1))
t = [0:0.01:5];
y = invgampdf(t/(n-1)/sample_var,(n-1)/2,1/2)/(n-1)/sample_var;
line([sample_var, sample_var],[0,max(y)])
%plot(t,y,'r')
kde(gibbs_sample(:,2)) % plot the histogram of the variance
%mean(gibbs_sample(:,2))
%sample_var
