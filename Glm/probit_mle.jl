using DelimitedFiles, Distributions, LinearAlgebra
affair = readdlm("affair.csv", ',')
y = affair[:,1];
X = affair[:,2:end];
n, k = size(X);
#find the MLE and the information matrix
S = ones(k,1); #score
betat = (X'*X)\(X'*y); #initial guess
e = 10^(-5);  #tolerance level
while sum(abs.(S)) > e  #stopping criterion
    Xbetat = X*betat;
    phi = pdf.(Normal(0,1),Xbetat);
    Phi = cdf.(Normal(0,1),Xbetat);
    global S = sum(repeat(y.*phi./Phi-(1 .- y).*phi./(1 .-Phi),outer = [1,k]).*X,dims=1)';
    d = phi.^2 ./ (Phi.*(1 .-Phi));
    IM = X'*diagm(vec(d))*X; # information matrix
    global betat = betat + IM\S;
end
println(round.(betat,digits=4))