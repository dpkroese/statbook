# CPI_ARMA.r
source("loglike_ARMA11.r");
##############<load a file>##############
file = as.matrix(read.csv("USCPI.csv"));#
file = matrix(file[,1]);                #
file = matrix(file[-258]);              #
d = 5.673853997; #The first line        #
USCPI = matrix(rbind(d,file),,1);       #
#########################################
y0 = USCPI[1];
y = matrix(USCPI[2:nrow(USCPI)]);
T = length(y);
X = cbind(matrix(rep(1,T)),rbind(y0,matrix(y[1:(nrow(y)-1)])));
f<-function(psi){
  res<-loglike_ARMA11(psi,X,y);
  result=-1*res[[1]];
  return (result);
}
psihat=optim(par=0, fn=f, gr = NULL, method = c("L-BFGS-B"),lower = -Inf, upper = Inf,control = list(), hessian = T)$par;
result <- loglike_ARMA11(psihat,X,y);
l <- result[[1]];
rhohat <- result[[2]];
sig2hat <- result[[3]];
print(l);
print(rhohat);
print(sig2hat);