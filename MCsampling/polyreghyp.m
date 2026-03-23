clear all
clf
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,3.8,4.8,1.7,...
    -0.4,4.5,1.3,0.4,2.6,4,2.9,1.6]';
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,...
    8.05,3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04]';
n = length(x);
X = [ones(n,1) x x.^2 x.^3 x.^4 x.^5];
m = size(X,2);
XX = X'*X;
Xy = X'*y;
betahat = XX\Xy; %estimate under the full model
ym = X*betahat; 

hold on %plotting data & 5th degree pol fit
plot(x,y,'.')
xvals = [min(x):0.1:max(x)]'; 
nn = size(xvals,1)
Xmat =  [ones(nn,1) xvals xvals.^2  xvals.^3 xvals.^4 xvals.^5];  
yvals =Xmat*betahat;
plot(xvals,yvals,'k')

X_4 = [ ones(n,1) x x.^2 x.^3 x.^4 ]; %omitting 5th degree
k = size(X_4, 2); %number of parameters in reduced model
betahat_4 = X_4'*X_4\(X_4'*y);
y_4 = X_4*betahat_4; 
T_4 = (n -m)/(m - k)*(norm(y-y_4)^2 - norm(y - ym)^2)/norm(y - ym)^2
pval_4 = 1 - cdf('F',T_4,m-k,n-m)



X = [ones(n,1) x x.^2 x.^3 x.^4 ];
m = size(X,2);
XX = X'*X;
Xy = X'*y;
betahat = XX\Xy; %estimate under the full model
ym = X*betahat; 
X_3 = [ ones(n,1) x x.^2 x.^3  ]; %omitting 4th degree
k = size(X_3, 2); %number of parameters in reduced model
betahat_3 = X_3'*X_3\(X_3'*y);
y_3 = X_3*betahat_3; 
T_3 = (n -m)/(m - k)*(norm(y-y_3)^2 - norm(y - ym)^2)/norm(y - ym)^2
pval_3 = 1 - cdf('F',T_3,m-k,n-m)


X = [ones(n,1) x x.^2 x.^3  ];
m = size(X,2);
XX = X'*X;
Xy = X'*y;
betahat = XX\Xy; %estimate under the full model
ym = X*betahat; 
X_2 = [ ones(n,1) x x.^2  ]; %omitting 3rd degree
k = size(X_2, 2); %number of parameters in reduced model
betahat_2 = X_2'*X_2\(X_2'*y);
y_2 = X_2*betahat_2; 
T_2 = (n -m)/(m - k)*(norm(y-y_2)^2 - norm(y - ym)^2)/norm(y - ym)^2
pval_2 = 1 - cdf('F',T_2,m-k,n-m)


Xmat =  [ones(nn,1) xvals xvals.^2  xvals.^3];  
yvals =Xmat*betahat;
plot(xvals,yvals,'r') %plotting 3rd degree pol fit
hold off

