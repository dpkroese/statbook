clear all
figure(1)
x = (rand(1,100) < 0.5); 
bar(x)
figure(2)
y = cumsum(x)./[1:100];
plot(y)
