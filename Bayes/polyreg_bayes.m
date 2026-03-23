% polyreg_bayes.m
x = [4.7,2,2.7,0.1,4.7,3.7,2,3.4,1.3,3.8,4.8,1.7,...
    -0.4,4.5,1.3,0.4,2.6,4,2.9,1.6]';
y = [6.57,5.15,7.15,0.18,6.48,8.95,5.24,10.54,1.24,...
    8.05,3.56,3.4,2.18,7.16,2.32,-0.23,7.68,9.09,9.13,4.04]';
n = length(x);
X = [ones(n,1) x x.^2 x.^3 x.^4 x.^5];
XX = X'*X;
Xy = X'*y;
m = 6;
N = 10^5;            % Gibbs sample size
V0 = eye(m)*100;     % prior for beta
invV0 = V0\eye(m);
alp0 = 2; lam0 = 1;  % prior for sig2
beta = XX\Xy;
sig2 = sum((y-X*beta).^2)/n;
gibbs_sample = zeros(N,m+1);
lpostden_sample = zeros(N,4); 
for k=1:N    
    D = (invV0 + XX/sig2)\eye(m);
    betahat = D*(Xy/sig2);
    beta = betahat + chol(D,'lower')*randn(m,1);
    sig2 = 1/gamrand(alp0+n/2,lam0+sum((y-X*beta).^2)/2);
    gibbs_sample(k,:)=[beta' sig2];
    lp1 = lmvnpdf([0 0 0 0]',betahat(3:end),D(3:end,3:end));
    lp2 = lmvnpdf([0 0 0]',betahat(4:end),D(4:end,4:end));
    lp3 = lmvnpdf([0 0]',betahat(5:end),D(5:end,5:end));
    lp4 = lmvnpdf(0,betahat(6),D(6,6));    
    lpostden_sample(k,:) = [lp1 lp2 lp3 lp4];
end
lpostden = zeros(4,1);
for k=1:4
    maxpden = max(lpostden_sample(:,k));
    lpostden(k) = log(mean(exp(lpostden_sample(:,k)-maxpden))) + maxpden;    
end
lpriden = zeros(4,1);
lpriden(1) = lmvnpdf([0 0 0 0]',[0 0 0 0]',V0(3:end,3:end));
lpriden(2) = lmvnpdf([0 0 0]',[0 0 0]',V0(4:end,4:end));
lpriden(3) = lmvnpdf([0 0]',[0 0]',V0(5:end,5:end));
lpriden(4) = lmvnpdf(0,0,V0(6,6));
lBF = lpostden - lpriden;

