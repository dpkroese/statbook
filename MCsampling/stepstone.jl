using Plots
N = 101; # times
P = [0    0.2  0   0.3  0.5  0;
     0.5  0    0.5 0    0    0;
     0.3  0    0   0.7  0    0;
     0    0    0   0    0    1;
     0    0    0   0.8  0    0.2;
     0    0    0.1 0    0.9  0];
x = zeros(Int64,N); x[1]= 1;
tot = zeros(6); tot[1] = 1;
for t in 1:N-1 # generate the Markov chain
    x[t+1] = minimum(findall(cumsum(P[x[t],:]) .> rand()));
    tot[x[t+1]] = tot[x[t+1]] + 1;
end
p = plot(0:N-1,x) # plot the path
scatter!(0:N-1,x)
println(tot/N) # fractions of visits to the states

using LinearAlgebra
f = nullspace(I - P');
f = f/sum(f)
