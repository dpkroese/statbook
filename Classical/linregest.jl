using LinearAlgebra,  Plots
N = 100; x = collect(1:N)/N ;
beta = [6; 13]; sigma = 2;# parameters
X = [ones(N,1) x];# design matrix
y = X*beta + sigma*randn(N); #draw the y-data
scatter(x,y) # plot the data
betahat = X'*X\(X'*y) # solve the normal equations
sigmahat = norm(y - X*betahat)/sqrt(N) # estimate for sigma