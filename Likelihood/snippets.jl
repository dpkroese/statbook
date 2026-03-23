using Plots
n = 16;
x = [5, 5, 4,4,4,4, 3,3,3,3,2,2,2,2,1,1];
sm = sum(x);
pr = prod(factorial.(x));
L(t) =t^sm/pr*(exp(t) -1)^(-n);
t = 0.001:0.001:5;
p1 = plot(t,L.(t))
p2 = plot(t,log.(L.(t)))
plot(p1,p2,layout=(2,1))
c,i = findmax(L.(t))
t[i]



using Roots
x = [1.1668, 0.0738, 0.7740,1.0160, 0.4822, 1.4559, 0.1752, 0.5209,
0.1537, 0.2947]
f(b) = sum(b./(b .- x)) .- 20;
bhat = fzero(f,1.7)

using StatsBase, SpecialFunctions
x = [
29.7679,   12.8406,  105.3225 ,  46.6101,   75.7135,   72.0340,
64.1004,   33.9008,   35.2510,   50.9201,   29.8086,   32.6963,
131.5229,  65.3381,   29.1369 ,  61.8774,   31.0650,   54.4877,
103.6889,  68.0230 ,  89.6879 ,  30.1994 ,  48.3140 ,  54.4447,
29.2253,   27.0242,  102.5929 ,  63.7344,   43.0354 ,  96.5552];
n = length(x);
sumlogx = sum(log.(x)); sumx = sum(x);
alpt = mean(x)^2/var(x); lamt = mean(x)/var(x); #intl. guess
thetat = [alpt, lamt]
for i=1:5 #just repeat the NR step 5 times
   global lamt, alpt, thetat
    S = [ n*(log(lamt) - digamma(alpt)) + sumlogx; n*alpt/lamt - sumx ];
    I = n * [trigamma(alpt) -1/lamt; -1/lamt alpt/lamt^2 ];
    thetat = thetat + inv(I)*S # using inv is OK (dim =2)
    alpt = thetat[1]; lamt = thetat[2];
end
print(thetat)


# using Newton-Raphson
S(t) =  34/t + 125/(t+2) - 38/(1-t);
H(t) = -34/t^2 - 125/(t+2)^2 - 38/(1-t)^2;
told = Inf;
t = 0.5; #initial guess
while abs(t - told)> 10^(-10)
told = t;
t = t - S(t)/H(t)
end


#using grid search
logL(t) = 125*log(2 + t) + 38*log(1-t) + 34*log(t);
t = 0.01:0.00001:0.99;
maxt, i = findmax(logL.(t));
t[i]


#Q6.23
using Plots
n=10;
a = 1.96;
T1 = x -> (a^2 + 2*x - a*sqrt(a^2-4*n*((x/n)-1).*(x/n)) )/
(2*(a^2 +n));
T2 = x ->  (a^2 + 2*x + a*sqrt(a^2-4*n*((x/n)-1).*(x/n)) )/
(2*(a^2 +n));
p = 0:0.01:1;
cvp = zeros(size(p));
for i=1:length(p);
tot = 0;
for x = 0:n
tot = tot + (T1(x) <= p[i] <= T2(x))*
binomial(n,x)*p[i]^x*(1-p[i])^(n-x);
end
cvp[i] = tot;
end
plot(p,cvp)