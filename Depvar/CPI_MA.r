# CPI_MA.r
source("loglike_MA1.r");
##############<load a file>##############
file = as.matrix(read.csv("USCPI.csv"));#
file = matrix(file[,1]);                #
file = matrix(file[-258]);              #
d = 5.673853997; #The first line        #
USCPI = matrix(rbind(d,file),,1);       #
#########################################
y0 = USCPI[1];
y = matrix(USCPI[2:nrow(USCPI)]);
Dely = y - rbind(y0,matrix(y[1:(nrow(y)-1)])); # define the dependent variable
T = length(Dely); # negative of the log-likelihood
f<-function(theta){
  res<-loglike_MA1(theta,Dely);
  return (-res);
}
theta0 = cbind(0, var(Dely));                # initial guess
thetahat=optim(par=theta0, fn=f, gr = NULL, method = c("L-BFGS-B"),lower = -Inf, upper = Inf,control = list(), hessian = T)$par;
psihat = thetahat[1];
l = loglike_MA1(thetahat,Dely);                
print(thetahat);
print(psihat);
print(l);