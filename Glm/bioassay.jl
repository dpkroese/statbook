using LinearAlgebra
y = [0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 1 1 1 1 1]';
x = repeat([-0.863 -0.296 -0.053 0.727], inner = (1,5));
X = [ones(20, 1) x'];        # design matrix
betat = (X' * X) \ (X' * y); # initial guess
S = ones(2, 1);              # score
IM = zeros(2,2)              # info matrix
e = 10^(-5); # tolerance level
while sum(abs.(S)) > e       # stopping criterion
    global betat, S, IM
    mu = 1 ./ (1 .+ exp.(-X * betat))
    S = sum(repeat((y - mu), outer=(1,2)) .* X, dims=1)'
    IM = X' * diagm(vec(mu .* (1 .- mu))) * X 
    betat = betat + IM \ S
end
V = IM \ I
println(betat)
println(V)