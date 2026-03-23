clear all
yy = [9.2988, 9.4978, 9.7604, 10.1025;
      8.2111, 8.3387, 8.5018,  8.1942;
      9.0688, 9.1284, 9.3484,  9.5086;
      8.2552, 7.8999, 8.4859,  8.9485];
n = numel(yy); [nrow,ncol] = size(yy);
y = yy(:);
X_1 = ones(n,1);
KM = kron(eye(ncol),ones(nrow,1)); X_2 = KM(:,1:ncol-1);
X_2(n-nrow+1:n,:) = -ones(nrow,ncol-1);
C = [eye(nrow-1), -ones(nrow-1,1)]'; X_3 = repmat(C,ncol,1);
X = [X_1 X_2 X_3];
m = size(X,2); %number of parameters in full model
betahat = X'*X\(X'*y); %estimate under the full model
ym = X*betahat; 
X_12 = [X_1 X_2]; %omitting the block effect
k = size(X_12, 2); %number of parameters in reduced model
betahat_12 = X_12'*X_12\(X_12'*y);
y_12 = X_12*betahat_12; 
T_12 = (n -m)/(m - k)*(norm(y-y_12)^2 - norm(y - ym)^2)/norm(y - ym)^2
pval_12 = 1 - cdf('F',T_12,m-k,n-m)
X_13 = [X_1 X_3]; %omitting the treatment effect
k = size(X_13, 2); %number of parameters in reduced model
betahat_13 = X_13'*X_13\(X_13'*y); 
y_13 = X_13*betahat_13;
T_13 = (n-m)/(m-k)*(norm(y-y_13)^2 - norm(y - ym)^2)/norm(y - ym)^2
pval_13 = 1 - cdf('F',T_13,m-k,n-m)

