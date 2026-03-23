%mcint.m
clear all
c = (2*pi)^(3/2);
H = @(z) c*sqrt(abs(sum(z,2))); 
N = 10^6; 
Z = randn(N,3); X = H(Z);
mX = mean(X); sX = std(X); 
R = 1.96*sX/sqrt(N);
fprintf('Estimate = %g, CI = (%g, %g)\n', mX, mX - R, mX + R)