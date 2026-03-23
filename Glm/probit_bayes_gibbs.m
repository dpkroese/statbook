% probit_bayes_gibbs.m
load('affair.csv');
y = affair(:,1); 
X = affair(:,2:end);
XX = X'*X;
[n k] = size(X);
V0 = 10*eye(k);            % prior covariance
invV0 = V0\speye(k);
burnin = 500;
nloop = 5000+burnin;
store_beta = zeros(nloop,k);
z = y;                     % initial guess
beta = XX\(X'*z);
                           % compute a few things before the loop
id0 = find(y==0); id1 = find(y==1);
n0 = length(id0); n1=n-n0;
V = (invV0 + XX)\speye(k); % posterior covariance
for i = 1:nloop
        % sample y^*
    Xb = X*beta;
    z(id0) = tnormrnd(Xb(id0),ones(n0,1),-inf,0);  
    z(id1) = tnormrnd(Xb(id1),ones(n1,1),0,inf); 
        % sample beta
    dbeta = X'*z;
    beta = V*dbeta + chol(V,'lower')*randn(k,1);
	store_beta(i,:) = beta';
end
store_beta = store_beta(burnin+1:end,:);  % discard the burnin
mean(store_beta)
std(store_beta)