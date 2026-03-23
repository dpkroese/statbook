% CPI_MA.m
load USCPI.csv
y0 = USCPI(1);
y = USCPI(2:end);
Dely = y - [y0; y(1:end-1)];           % define the dependent variable
T = length(Dely);
                                       % negative of the log-likelihood
f = @(theta) -loglike_MA1(theta,Dely);
theta0 = [0 var(Dely)];                % initial guess
thetahat = fminsearch(f, theta0)
psihat = thetahat(1)
l = loglike_MA1(thetahat,Dely)