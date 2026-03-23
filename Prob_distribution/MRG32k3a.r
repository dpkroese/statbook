m1=2^32-209; m2=2^32-22853;
ax2p=1403580; ax3n=810728;
ay1p=527612; ay3n=1370589;

X=array(12345,3); 
Y=array(12345,3);

N=100;
U=matrix(1,1,N);
for (t in 1:N)
  {
    Xt=(ax2p*X[2]-ax3n*X[3])%%m1;
    Yt=(ay1p*Y[1]-ay3n*Y[3])%%m2;
    { 
        if(Xt<=Yt)
            {
              U[t]=(Xt-Yt+m1)/(m1+1);
             } 
             else
               {
               U[t]=(Xt-Yt)/(m1+1);
               }
    }
    X[2:3]=X[1:2]; X[1]=Xt; Y[2:3]=Y[1:2];Y[1]=Yt;
}