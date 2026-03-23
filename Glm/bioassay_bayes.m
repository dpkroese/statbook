% bioassay_bayes.m

y = [0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 1 1 1 1 1]';
x = repmat([-0.863 -0.296 -0.053 0.727],5,1);
A = [ones(20,1) x(:)];      % design matrix
betat = (A'*A)\(A'*y);      % initial guess
S = ones(2,1);              % score
e = 10^(-5);                % tolerance level
while sum(abs(S)) > e       % stopping criteria
    mu = 1./(1+exp(-A*betat));
    S = sum(repmat((y - mu),1,2).*A)';  
    I = A'*diag(mu.*(1-mu))*A; %info matrix 
    betat = betat + I\S;  
end
V = I\eye(2);

B = chol(V,'lower');
burnin = 500;
nloop = 10000+burnin; 
store_beta = zeros(nloop,2);
nu = 5;			     % df for the proposal 
% log posterior density
logf = @(b) sum((y-1).*(A*b)-log((1+exp(-A*b))));
% log density of the t proposal
logprop = @(b) -.5*(nu+2)*log(1+(b-betat)'*(V\(b-betat))/nu);
beta = betat;               % initialize the chain
for i = 1:nloop
			     % candidate draw from the t proposal
    betac = betat + B*randn(2,1)*sqrt(nu/gamrand(nu/2,1/2)); 
    rho = logf(betac)-logf(beta) + ...
        +logprop(beta)-logprop(betac);
    if exp(rho) > rand
        beta = betac;
    end
    store_beta(i,:) = beta';    
end
store_beta = store_beta(burnin+1:end,:);  % discard the burnin
mean(store_beta)
cov(store_beta)
