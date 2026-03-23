%genlink.m
x = [ 125 18 20 34 ]';
n = sum(x); 
theta = 4*(x(1)/n-1/2);    %% initial guess
err = 1;                   
while abs(err) > 10^(-5)    %% stopping criteria
    z = 2*x(1)/(2+theta);         %% E-step
    temp = (x(1)+x(4) - z)/(n-z); %% M-step
    err = theta - temp;       
    theta = temp    
end

