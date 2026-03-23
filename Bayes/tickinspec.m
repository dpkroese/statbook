%tickinspec.m
clear all,clf
c = 3.481048347*10^(19);
t = 0:0.01:4;
f = t.^(60).*(3^(59)*exp(-92*t/3) + 2^59*exp(-21*t) + exp(-12*t))/c;
hold on
p = plot(t,f,'k','LineWidth',2);
% Uncomment the following line to preserve the X-limits of the axes
%  xlim([1 4]);
% % Uncomment the following line to preserve the Y-limits of the axes
%  ylim([0 1.2]);
% hold('all');
% 
% % Create plot
% plot(t,f,'LineWidth',2,'Color',[0 0 0]);
% 
% % Create xlabel
% xlabel('$t$','Interpreter','latex','FontSize',20);
% 
% % Create ylabel
% ylabel('$f(t \, |\, x = 60)$','Interpreter','latex','FontSize',20);
%     

n = 10000;
x = 60; %number of tickets
p = [1/3,1/3,1/3] %initial value
kk = zeros(1,n);
tt = zeros(1,n);
k = min(find(cumsum(p)> rand));
for i=1:n
    a = x+ 1;
    b = 2/k + 10*k; 
    t = gamrand(a,b);
    tt(i) = t;
    p = [exp(-12*t),2^(x-1)*exp(-21*t), 3^(x-1)*exp(-92*t/3)];
    p = p/sum(p);
    k = min(find(cumsum(p)> rand));
    kk(i) = k;
end
p1est = sum(kk == 1)/n  %estimate of posterior probability 1
p2est = sum(kk == 2)/n  %estimate of posterior probability 2
p3est = sum(kk == 3)/n  %estimate of posterior probability 3
kde(tt)
