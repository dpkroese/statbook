% probit_mle.m

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
    I = X'*sparse(1:n,1:n,d)*X; % infomation matrix
    betat = betat + I\S;    
end
betat