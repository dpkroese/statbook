% probit_bayes.m
load('affair.csv');
y = affair(:,1);
X = affair(:,2:end);
[n k] = size(X);    
    %% obtain the proposal density for the MH step
S = ones(2,1);              % score
betat = (X'*X)\(X'*y);      % initial guess
e = 10^(-5);                % tolerance level
while sum(abs(S)) > e       % stopping criteria
    Xbetat = X*betat;
    phi = normpdf(Xbetat);  Phi = normcdf(Xbetat);
    S = sum(repmat(y.*phi./Phi-(1-y).*phi./(1-Phi),1,k).*X)'; 
    d = phi.^2./(Phi.*(1-Phi));
    I = X'*sparse(1:n,1:n,d)*X;
    betat = betat + I\S;    
end
burnin = 500;
nloop = 5000+burnin;
V = I\eye(k);               % scale matrix for the proposal
B = chol(V,'lower'); 
nu = 5;                     % df for the proposal
b0 = zeros(k,1);            % prior mean
V0 = 10*eye(k);             % prior covariance
% log-posterior density
logf=@(b) y'*log(normcdf(X*b)) + (1-y)'*log(normcdf(-X*b)) ...
    - .5*(b-b0)'*(V0\(b-b0));
% log-proposal density
logprop=@(b) -.5*(k+nu)*log(1+(b-betat)'*(V\(b-betat))/nu);
store_beta = zeros(nloop,k);
beta = betat;
for i = 1:nloop
        % candidate draw from the t proposal
    betac=betat + B*randn(k,1)*sqrt(nu/gamrand(nu/2,1/2));
    rho = logf(betac)-logf(beta) + ...
        + logprop(beta)-logprop(betac);
    if exp(rho) > rand
        beta = betac;
    end
    store_beta(i,:) = beta';
end
store_beta = store_beta(burnin+1:end,:); % discard the burnin
mean(store_beta)