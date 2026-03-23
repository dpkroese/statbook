using Plots
f(x,y) = exp(-4*(y-x^2)^2 + (y-1)^2)*(y < 2)
N = 10000
xx = zeros(N,2); x = [0 -1]; xx[1,:] = x;
for i in 2:N
    y = x + randn(1,2); #proposal
    alpha = min(f(y[1],y[2])/f(x[1],x[2]),1); #acceptance prob
    r = (rand() < alpha);
    global x = r*y + (1-r)*x; #next value of the Markov chain
    xx[i,:] = x;
end
scatter(xx[:,1],xx[:,2],markersize = 1)
x = range(-2, stop=2, length=50)
y = range(-2, stop=2, length=50)
contour!(x,y,f)