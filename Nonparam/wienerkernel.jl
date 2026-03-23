using LinearAlgebra, Plots
x = 0:1:2*pi
n = length(x)
y = sin.(x)  
k(x,u) = min(x,u) # kernel
#k(x,u) = sinc(x-u)  #sinc kernel
K = zeros(n,n)
for i=1:n
    for j=1:n
        K[i,j] = k(x[i], x[j]) 
    end
end
alpha = pinv(K)*y     # compute an optimal alpha
xx = 0:0.01:2*pi
N = length(xx); g = zeros(N);
Kx = zeros(n,N)
for i=1:n
    for j=1:N
       Kx[i,j] = k(x[i],xx[j]) 
    end
end
g =  Kx'*alpha;  # function values

plot(xx,sin.(xx),color=:black,linestyle=:dash)
plot!(xx,g,color=:blue)
scatter!(x,y,legend=false,color=:black)
