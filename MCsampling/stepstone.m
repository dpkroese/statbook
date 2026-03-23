%stepstone.m
clear all, clf
n = 101;
P = [0,   0.2, 0,   0.3, 0.5, 0;
     0.5, 0,   0.5, 0,   0,   0;
     0.3, 0,   0,   0.7, 0,   0;
     0,   0,   0,   0,   0,   1 ;
     0,   0,   0,   0.8, 0,   0.2;
     0,   0,   0.1, 0,   0.9, 0];
x = zeros(1,n); x(1)= 1;
tot = zeros(1,6); tot(1) = 1;
for t=1:n-1                            % generate the Markov chain
    x(t+1) = min(find(cumsum(P(x(t),:))> rand));
    tot(x(t+1)) = tot(x(t+1)) + 1;
end
hold on,  plot(0:n-1,x,'.'), plot(0:n-1,x),  hold off  % plot the path
tot/n
f = null(eye(6) - P')';
f = f/sum(f)