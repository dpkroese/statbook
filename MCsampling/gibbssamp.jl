using Plots
f(x,y) = exp(-(x*y + x + y))*(x > 0 && y > 0)
N = 10^4; x = zeros(N,2); x2 = 1;
for i  in  2:N
    x1 = -log(rand())/(x2+1);
    global x2 = -log(rand())/(x1+1);
    global x[i,:] = [x1 x2];
end
scatter(x[:,1],x[:,2],markersize=0.2)