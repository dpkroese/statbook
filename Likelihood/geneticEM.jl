x = [ 125 18 20 34 ]; n = sum(x);
theta = 4*(x[1]/n-1/2);
# initial guess
err = 1;
while abs(err) > 10^(-5)
    # stopping criterion
    z = 2*x[1]/(2+theta);
    # E-step
    temp = (x[1]+x[4] - z)/(n-z); # M-step
    global err = theta - temp;
    global theta = temp;
end