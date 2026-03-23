%linregest.m
clear all
N = 100;
x = ( (1:N)/N )';
beta = [6;13];
sigma = 2;
A = [ones(N,1),x];
y = A*beta + sigma*randn(N,1);
plot(x,y,'.');
betahat = A'*A\(A'*y)
%betahat2 = pinv(A)*y  %the same
sigmahat = norm(y - A*betahat)/sqrt(N)
