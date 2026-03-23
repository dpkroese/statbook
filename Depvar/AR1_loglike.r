# AR1_loglike.r
AR1_loglike<-function(rho,y,X) {
T = length(y);
#### Making a sparse matrix#####
Sp=diag(T);                    #
d=matrix(rep(0,T),1);          #
sparse=rbind(d,Sp);            #
sparse<-sparse[-(T+1),]        #
################################
H = diag(T) - rho*sparse;
HH = t(H) %*% H;
betahat = solve((t(X) %*% HH %*% X),(t(X) %*% HH %*% y));
e = y-X %*% betahat;
sigma2hat = t(e) %*% HH %*% e/T;
l = -T/2*log(2*pi*sigma2hat) - .5/sigma2hat*t(e) %*% HH %*% e;
result=list(l,betahat,sigma2hat);
return (result);
}