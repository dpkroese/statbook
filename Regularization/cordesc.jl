using LinearAlgebra,StatsBase, Plots,LaTeXStrings
beta = [0.1, 0.2, 0.3]
X = [8.46  4.73   6.26;
     5.29  0.98  5.54;
     6.99  2.96  0.38;
     9.87  4.67  8.89;
     9.58  9.22  5.98]
n,p = size(X)
r = [-0.83, 0.93, -0.24, -0.40, -0.02]
y = X*beta + r;

Xmu=mean(X,dims=1); ymu=mean(y);
y=y .- ymu; X=X .- Xmu;  # centering 

function ADMM(X,y,lam)
   rho  = 100.0; eps = 10e-9
   kappa = lam / (2 * rho)
   betahat = (X' * X) \ (X' * y)
   b = z = u = betahat
   cont = true
   while cont
      bold = b
      b = (X' * X + rho * I) \ (X' * y + rho * (z - u))
      a = b + u
      z = max.(a .- kappa, 0) - max.(-a .- kappa, 0)
      u = u + b - z
      cont = norm(b - bold) > eps
   end
   return b
end

lams = 0:0.1:30
bs = zeros(p, length(lams))
norms = zeros(length(lams))
for i = 1:length(lams)
   lambda = lams[i]
   b = ADMM(X,y,lambda)
   bs[:, i] = b
   norms[i] = norm(b, 1)
end
plot(tickfontsize=15,guidefontsize=20,tickfont = "Computer Modern")
xlabel!(L"\mathrm{1-norm}")
ylabel!(L"\bf{\beta}")
for pp=1:p
   global g
   if beta[pp] != 0
     g = plot!(norms,bs[pp,:],legend=false,color="blue",linewidth=2)
   else
     g = plot!(norms,bs[pp,:],legend=false,color="orange",linestyle=:dash,linewidth=2)
   end
   display(g)
end
g = plot!([0,1.3],[0,0],color=:black)

savefig(g,"admmex.pdf")