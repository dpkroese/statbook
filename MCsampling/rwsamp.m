%rwsamp.m
clear all, clc, clf
f=@(X,Y)exp(-4*(Y-X.^2).^2+(Y-1).^2).*(Y < 2);
[X,Y] = meshgrid(-2:0.05:2,-2:0.05:2);
Z = f(X,Y);
figure(1)
surfl(X,Y,Z)
shading interp
colormap gray
N = 10000; %sample size
xx = zeros(N,2); x = [0,-1]; xx(1,:) = x;
for i=2:N 
    z =  randn(1,2); %do standard gaussian first
    y = x + z;
    alpha = min(f(y(1),y(2))/f(x(1),x(2)),1); %acceptance prob
    r = (rand < alpha);
    x = r*y + (1-r)*x; %new x-value
    xx(i,:) = x;
    
end
figure(2)
hold on
plot(xx(:,1),xx(:,2),'ko','MarkerSize',1)
contour(X,Y,Z,'LineColor','k','LineWidth',2)