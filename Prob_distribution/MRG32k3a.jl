m1 = 2^32 - 209
m2 = 2^32 - 22853
ax2p = 1403580
ax3n = 810728
ay1p = 527612
ay3n = 1370589

X = [12345, 12345, 12345]  # Initial X at -1, -2, -3
Y = [12345, 12345, 12345]  # Initial Y at -1, -2, -3

N = 1000  # Compute the sequence for N steps
U = zeros(Float64, N)

for t in 1:N
    global X, Y
    Xt = mod(ax2p * X[2] - ax3n * X[3], m1)
    Yt = mod(ay1p * Y[1] - ay3n * Y[3], m2)
    
    if Xt <= Yt
        U[t] = (Xt - Yt + m1) / (m1 + 1)
    else
        U[t] = (Xt - Yt) / (m1 + 1)
    end
    
    X = [Xt, X[1], X[2]]
    Y = [Yt, Y[1], Y[2]]
end
