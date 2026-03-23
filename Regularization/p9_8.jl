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


S(x,gam) = x*max(1.0 - gam/abs(x),0)

function CD(X, y, lam)
   z = diag(X' * X)
   beta = (X' * X) \ (X' * y)
   u = y - X*beta
   for loop = 1:10000  
      bold = copy(beta)  # important
      for j = 1:p
         b = S(beta[j]+ u'*X[:,j]/z[j], lam/(2*z[j]))
         u = u + (beta[j]-b)*X[:,j]
         beta[j] = b  
      end
      if norm(beta-bold)<10e-12
         break
      end
   end
   return beta
end
lam = 10
print(CD(X,y,lam))

lams = 0:0.1:30
bs = zeros(p, length(lams))
norms = zeros(length(lams))
for i = 1:length(lams)
   lambda = lams[i]
   b = CD(X,y,lambda)
   bs[:, i] = b
   norms[i] = norm(b, 1)
end
plot(tickfontsize=15,guidefontsize=20,tickfont = "Computer Modern")
xlabel!(L"\mathrm{1-norm}")
ylabel!(L"\bf{\beta}")
for pp=1:p
      display(plot!(norms,bs[pp,:],legend=false,color="blue",linewidth=2))
end
plot!([0,1.3],[0,0],color=:black)