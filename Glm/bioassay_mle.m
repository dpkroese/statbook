% bioassay.m
y = [0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 1 1 1 1 1]';
x = repmat([-.863 -.296 -.053 .727],5,1);
X = [ones(20,1) x(:)];      % design matrix
betat = (X'*X)\(X'*y);      % initial guess
S = ones(2,1);              % score
e = 10^(-5);                % tolerance level
while sum(abs(S)) > 10^(-5) % stopping criteria
    mu = 1./(1+exp(-X*betat));
    S = sum(repmat((y - mu),1,2).*X)';
    I =   X'*sparse(1:20,1:20,mu.*(1-mu))*X; % infomation matrix
    betat = betat + I\S;    
end
V = I\eye(2);               % covariance matrix
