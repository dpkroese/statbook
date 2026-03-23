using SparseArrays, LinearAlgebra
y = [8 16 15 23 5 13 4 19 33 2]'
birth = [236 739 970 2371 309 679 26 1272 3246 28]'
Htype = [0 1 1 1 1 1 0 1 1 0]';
n = length(y)
X = [ones(n,1) birth Htype];
betat = (X'*X)\(X'*log.(y .+ .001)) # initial guess
S = ones(2,1)                       # score
e = 10^(-5)                         # tolerance level
IM = zeros(2,2)
while sum(abs.(S)) > 10^(-5)        # stopping criterion
    global S, betat, IM
    mu = exp.(X*betat)
    S = sum(repeat((y - mu), 1,3).*X,dims=1)'  
    IM =   X'*sparse(1:n,1:n,vec(mu))*X # info matrix
    betat = betat + IM\S
end
println(round.(betat,digits=4))
V = IM\I  # inverse of the info matrix
println(round.(V,digits=10))