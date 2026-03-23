# sales.r
source("AR1_loglike.r");
############<load a file>#############
ads = as.matrix(read.csv("ads.csv"));#
d = c(23.45,5.41); #The first line   #
ads = matrix(rbind(d,ads),ncol=2);   #
######################################
y = matrix(ads[,1],,1);
T = length(y);
X = cbind(matrix(rep(1,T)), ads[,2]);
f<-function(rho){
  res<-AR1_loglike(rho,y,X);
  result=-1*res[[1]];
  return (result);
}
rhohat=optimize(f,c(0,100))$minimum;
result<-AR1_loglike(rhohat,y,X);
l<-result[[1]];
betahat<-result[[2]];
sigma2hat<-result[[3]];
e = y - X %*% betahat;
#### Making a sparse matrix#####
Sp=diag(T);                    #
d=matrix(rep(0,T),1);          #
sparse=rbind(d,Sp);            #
sparse<-sparse[-31,]           #
################################
H = diag(T) - rhohat*sparse;
u = H %*% e;
#win.graph() #------------>For Windows users
quartz()#---------------->For Mac users
plot(u,xlab='',ylab='');
lines(c(0,T),c(0,0),col='blue');