x = [1.3615 3.5616 -14.2411 -4.4950 2.3014 1.1066 -9.3409 0.3779 0.9386 -0.1838]; #the data
a = 2; #initial guess
n = 10;
for i=1:7 
  println(a)
   a = a + 4*sum( (x .- a)./(1 .+ (x .- a).^2) )/n # note vectorization!
end