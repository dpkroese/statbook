%linregestconf.m
clear all, clf
N = 100; x = ( (1:N)/N )';
beta = [6;13]; sigma = 2;
X = [ones(N,1),x];
y = X*beta + sigma*randn(N,1);
plot(x,y,'.');
betahat = X'*X\(X'*y)
sigmahat = norm(y - X*betahat)/sqrt(N)
tquant = icdf('t',0.975,N-2)
ucl = zeros(1,N); lcl = zeros(1,N);
rl = zeros(1,N); el = zeros(1,N);
u=0;
for i=1:N
    u = u + 1/N;
    a = [1;u];
    rl(i) = a'*beta;
    ucl(i) = a'*betahat + tquant*norm(y - X*betahat)* ...
        sqrt(a'*inv(X'*X)*a)/sqrt(N-2);
    lcl(i) = a'*betahat - tquant*norm(y - X*betahat)* ...
        sqrt(a'*inv(X'*X)*a)/sqrt(N-2);
end
hold on
plot(x,rl), plot(x,ucl),plot(x,lcl)
hold off
