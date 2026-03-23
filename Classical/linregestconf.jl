using Distributions, Random, Statistics, LinearAlgebra, Plots
N = 100; x =  collect(1:N)/N ;
beta = [6; 13]; sigma = 2;# parameters
X = [ones(N,1) x];# design matrix
y = X*beta + sigma*randn(N,1); #draw the y-data
scatter(x,y); # plot the data
betahat = X'*X\(X'*y)  # solve the normal equations
sigmahat = norm(y - X*betahat)/sqrt(N) #estimate for sigma

tquant = quantile(TDist(N-2),0.975) # 0.975 quantile
ucl = zeros(N); lcl = zeros(N); # upper/lower conf. limits
rl = zeros(N)  # (true) regression line
u=0
for i in 1:N
   global u = u + 1/N
   a = [1;u]
   rl[i] = a'*beta;
   ucl[i] = dot(a,betahat) .+ tquant*norm(y - X*betahat)*sqrt(a'*inv(X'*X)*a)/sqrt(N-2)
   lcl[i] = dot(a,betahat) .- tquant*norm(y - X*betahat)*sqrt(a'*inv(X'*X)*a)/sqrt(N-2)
end
 plot!(x,rl,legend=false); plot!(x,ucl,legend=false); plot!(x,lcl,legend=false); scatter!(x,y,legend=false)