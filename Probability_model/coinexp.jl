x = rand(100) .< 0.5
t = 1:100
using Plots
bar(t,x,legend=false,color=:darkblue)
y = cumsum(x)./t
plot(t,y)