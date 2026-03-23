rm(list=ls())
win.graph()
x=(runif(100))<0.5
barplot(x)
win.graph()
y=cumsum(x)/(1:100);
plot(y, type='l')
