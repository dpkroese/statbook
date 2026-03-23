%indepsamp.m
clear all, clf
N = 10^5;
xx = zeros(1,N);
lam = 1;
f = @(x) x.^2.*exp(-x.^2 + sin(x));
g = @(x) lam*exp(-abs(x)*lam)/2;
alpha = @(x,y) min(f(y)*g(x)/(f(x)*g(y)), 1);
x = 0; xx(1) = x;
for t = 2:N+1
    y = -log(rand)*(2*(rand < 1/2) - 1);
    if rand < alpha(x,y)
        x = y;
    end
    xx(t) = x;
end
hold on
kde(xx(1:N+1));
c = quad(f,-5,5);
tt = [-4:0.1:4]; plot(tt,f(tt)/c,'r')
hold off