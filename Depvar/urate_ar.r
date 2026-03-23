# urate_ar.r
urate = matrix(c(5.7,5.8,5.7,5.9,5.9,6.1,6.1,5.8,5.7,5.6,5.4,5.4,5.3,5.1,5.0,5.0,4.7,4.6,4.6,4.4,4.5,4.5,4.7,4.8, 5.0,5.3,6.0,6.9,8.3,9.3,9.6,9.9,9.8,9.6,9.5,9.6,9.0,9.0,9.1,8.7));
y = matrix(urate[3:nrow(urate)]);
T = length(y);
X = cbind(matrix(rep(1,T)), matrix(urate[2:(nrow(urate)-1)]));
rhohat = solve((t(X)%*%X),(t(X)%*%y));
yhat1 = X %*% rhohat;                 # fitted values
uhat = y-yhat1;                   # residuals
sig2hat = t(uhat) %*% uhat/T;
#win.graph() #------------>For Windows users
quartz()#---------------->For Mac users
plot(y,type='l',lty=2,col='blue',xlab='',ylab='');
lines(yhat1,col='blue');