%gibbssamp.m
clear all
clc, clf
f=@(X1,X2)exp(-(X1.*X2 + X1 +X2)).*(X1 > 0 & X2 > 0);
N = 10^4; %sample size
x = zeros(N,2); x2 = 1;
for i=2:N 
    x1 = -log(rand)/(x2 +1);
    x2 = -log(rand)/(x1+1);
    x(i,:) = [x1,x2];
end
plot(x(:,1),x(:,2),'.')
