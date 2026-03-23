% sfran2.m
y = [1.39 1.29 1.12 1.16 1.52 ...
     1.62 1.88 1.87 1.24 1.18 ...
     0.95 0.96 0.82 0.92 1.18 ...
     1.20 1.47 1.41 1.57 1.65]';
Xalpha = kron(speye(5),ones(4,1));
Xgamma = kron(speye(10),ones(2,1));
f = @(theta) -sfran2_loglike(theta(1),theta(2),theta(3), ...
    theta(4),y,Xalpha,Xgamma);
yhat = mean(reshape(y,4,5));
theta0 = [mean(y) log(var(yhat)) log(var(y)/3) log(var(y)/3)];
thetahat = fminsearch(f,theta0);
thetahat(1)
exp(thetahat(2:4))