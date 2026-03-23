using LinearAlgebra,StatsBase, Plots

p = 20
n = 100
beta = collect(0.1*(1:p)) 
beta[div(p,2):p] = zeros(div(p,2)+1);

sig=3
X = 10*rand(n,p);
y = X*beta + sig*randn(n);

Xmu = mean(X, dims=1)
ymu = mean(y)
y = y .- ymu
X = X .- Xmu
X = X ./ std(X,dims=1)

function ADMM(X,y,lam)
   rho  = 0.1
   kappa = lam / (2 * rho)
   betahat = (X' * X) \ (X' * y)
   b = betahat*1000
   z = betahat*1000
   u = betahat*1000
   cont = true
   while cont
      bold = b
      b = (X' * X + rho * I) \ (X' * y + rho * (z - u))
      a = b + u
      z = max.(a .- kappa, 0) - max.(-a .- kappa, 0)
      u = u + b - z
      cont = norm(b - bold) > 10e-9
      #println(norm(b  - bold))
   end
   return b
end

rho = 0.1
lams = 0:0.05:3;
bs = zeros(p, length(lams));
norms = zeros(length(lams));
for i = 1:length(lams)
   println(i)
   lambda = lams[i]
   b = ADMM(X,y,lambda)
   bs[:, i] = b
   norms[i] = norm(b, 1)
end
#for pp=1:p
   #if beta[pp]
   #%plot(lams,bs(pp,:),"Color",C(3,:),'LineWidth',2)
   pp=10
   plot!(norms,bs[pp,:])
   #else
    #%   plot(lams,bs(pp,:),"Color",C(1,:),'LineWidth',2)
    #plot(norms,bs(pp,:),"Color",C(1,:),'LineWidth',2)
   #end
#end