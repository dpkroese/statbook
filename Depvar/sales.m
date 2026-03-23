% sales.m
load 'ads.csv';
y = ads(:,1);
T = length(y);
X = [ones(T,1) ads(:,2)];
f = @(rho) -AR1_loglike(rho,y,X);
rhohat = fminsearch(f,0);
[l betahat sigma2hat] = AR1_loglike(rhohat,y,X);

e = y - X*betahat;
H = speye(T) - rhohat*sparse(2:T,1:T-1,ones(1,T-1),T,T);
u = H*e;
hold on
plot(u,'o','MarkerSize',10,'MarkerEdgeColor','k'); 
line([0,T],[0,0]);