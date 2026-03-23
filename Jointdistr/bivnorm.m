%bivnorm.m
N = 1000;
rho = 0.8;
Sigma = [1 rho; rho 1];
B = chol(Sigma,'lower');
x = B*randn(2,N);
plot(x(1,:),x(2,:),'.');