using StatsBase, Distributions
x = [1.04, 3.30, 1.40, 1.53, 1.01, 1.02, 3.93, 1.61, 1.93, 2.25]
y = [1.00, 1.40, 1.30, 3.95, 0.08, 1.33, 0.66, 0.73, 1.49, 0.43]
m = length(x); n = length(y); N = m+n;
global count = 0
for i=1:m
   for j = 1:n
      if x[i] > y[j]
        global count = count + 1;
      end
   end
end
println(count) 
println(131 - m*(m+1)/2)