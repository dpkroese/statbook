using Plots, LinearAlgebra, Random
N = 1000; rho = 0.8;
Sigma = [1 rho; rho 1];
B = cholesky(Sigma).L;
# lower-triangular Cholesky matrix
x = B*randn(2,N);
scatter(x[1,:],x[2,:],ms=2,msw=0,legend=false)