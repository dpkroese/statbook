function SVRW(ystar, h, omega2h, Vh)
    T = length(h)
    # parameters for the Gaussian mixture
    pi = [0.0073 0.10556 0.00002 0.04395 0.34001 0.24566 0.2575]
    mui = [-10.12999 -3.97281 -8.56686 2.77786 0.61942 1.79518 -1.08819]
    sig2i = [5.79596 2.61369 5.17950 0.16735 0.64009 0.34023 1.26261]
    sigi = sqrt.(sig2i)
    s = zeros(Int64, T)
    for t = 1:T
        q = zeros(7)
        for i = 1:7
            q[i] = pi[i]*pdf(Normal(mui[i]-1.2704+h[t],sigi[i]),ystar[t])
        end
        q = q ./ sum(q)
        s[t] = minimum(findall(cumsum(q) .> rand()))
    end
    H = sparse(I, T, T) - sparse(2:T,1:(T-1),ones(T - 1), T, T)
    invOmegah = spdiagm(vec([1/Vh; 1/omega2h*ones(T - 1, 1)]))
    d = mui[s]
    invSigystar = spdiagm(vec(1 ./ sig2i[s]))
    Kh = H' * invOmegah * H + invSigystar
    R = cholesky(Kh)
    P = sparse(1:T, R.p, ones(T))
    Ch = P' * sparse(R.L)
    hhat = Kh \ (invSigystar * (ystar - d))
    h = hhat + Ch' \ randn(T, 1)
    return h, s
end
