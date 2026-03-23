using SparseArrays, LinearAlgebra, StatsBase,Plots,DelimitedFiles
y = readdlm("USCPI.csv")
T = length(y)
omega2_0 = 9      # initial condition
omega = .5^2      # fix omega
H = sparse(I,T,T) - sparse(2:T,1:(T-1),ones(T-1),T,T)
invOmega = sparse(1:T,1:T,vcat(1/omega2_0, 1/omega*ones(T-1)),T,T)
HinvOmegaH = H'*invOmega*H
sigma2t = var(y)  # initial guess
err = 1
while err> 10^(-4)
        # E-step
    Kt = HinvOmegaH + sparse(I,T,T)/sigma2t
    taut = Kt\(y/sigma2t)
        # M-step
    lam =  eigvals(Matrix(Kt))
    newsigma2t = (sum(1 ./ lam) .+ only((y-taut)'*(y-taut)))/T
        # update
    err = abs(sigma2t-newsigma2t)
    sigma2t = newsigma2t
end
Kt = HinvOmegaH + sparse(I,T,T)/sigma2t
taut = Kt\(y/sigma2t)
plot(y)
plot!(taut)
