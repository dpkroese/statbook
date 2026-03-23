using SparseArrays, LinearAlgebra
y = [6 3 2 1 3 1 1 8 2 0 2 6 1 3]'
RD = [26 21 19 11 21 16 19 29 13 5 3 29 3 21]'
n = length(y)
X = [ones(n,1) RD]
betat = (X'*X)\(X'*log.(y .+ .001)) # initial guess
S = ones(2,1)               # score
e = 10^(-5)  # tolerance level
IM = zeros(2,2)
while sum(abs.(S)) > 10^(-5) # stopping criterion
    global S, betat, IM
    mu = exp.(X*betat)
    S = sum(repeat((y - mu), 1,2).*X,dims=1)'  
    IM =   X'*sparse(1:n,1:n,vec(mu))*X # info matrix
    betat = betat + IM\S
end
println(round.(betat,digits=4))
V = IM\I  # inverse of the info matrix
println(round.(V,digits=4))