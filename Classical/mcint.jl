using Statistics, Printf
c = (2*pi)^(3/2); N = 10^6;
H = z -> c*sqrt.(abs.(sum(z,dims=2)))
Z = randn(N,3); X = H(Z);
mX = mean(X); sX = std(X);
R = 1.96*sX/sqrt(N);
LCI = mX - R; UCI = mX + R;
@printf("Estimate = %.3f, CI = (%.3f,%.3f)",mX,LCI,UCI)