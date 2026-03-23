%confintnorm.m
clear all
mu = 3; sig = 0.5; % parameters
alpha = 0.05; n = 10;
mu_count = 0; sig_count = 0;
for k=1:100
    x = mu + randn(n,1)*sig; % draw the iid sample
    mu_est = mean(x);  % estimate mu
    sig_est = std(x);  % estimate sigma
    tq = icdf('t',1-alpha/2,n-1);
    mu_lo = mu_est - tq*sig_est/sqrt(n); %low bound CI for mu
    mu_hi = mu_est +  tq*sig_est/sqrt(n); %upper bound 
    cq1 = icdf('chi2',1-alpha/2,n-1);
    cq2 = icdf('chi2',alpha/2,n-1);
    sig_lo = (n-1)*sig_est^2/cq1; %lower bound CI for sigma
    sig_hi = (n-1)*sig_est^2/cq2;  % upper bound
    mu_count = mu_count + (mu > mu_lo & mu < mu_hi);
    sig_count =  sig_count + (sig^2 > sig_lo & sig^2 < sig_hi);
end
disp([mu_count, sig_count]) % final counts
