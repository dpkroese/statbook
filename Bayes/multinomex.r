#multinomex.m
rm(list=ls());
source("kde.r");
source("dirichrnd.r");
x = matrix(c(53,57,147,93,38,113),1);
N = 10000;
p = array(0,c(N,2,3));
a = array(0,c(N,2,3));
p_row=array(0,c(2,N));
p_col=array(0,c(3,N));
alpha = x + 1;
for (i in 1:N){
  r = matrix(dirichrnd(alpha),1);
  h = cbind(r,1-sum(r));
  p[i,,] = t(array(t(h),c(3,2)));
}
for (i in 1:2){
  p_row[i,] = matrix(apply(p[,i,],1,sum),1);
}
for (j in 1:3){
  p_col[j,] = matrix(apply(p[,,j],1,sum),1);
}
for (i in 1:2){
  for (j in 1:3){
    a[,i,j] = p[,i,j] - t(p_row[i,])*t(p_col[j,]);  
  }
}
for (j in 1:3){
  result<-kde(a[,1,j],len=0,loop=j);
}