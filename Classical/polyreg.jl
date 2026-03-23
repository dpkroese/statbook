using LinearAlgebra
x = [4.7 2 2.7 0.1 4.7 3.7 2 3.4 1.3 3.8 4.8 1.7 -0.4 4.5 1.3 0.4 2.6 4 2.9 1.6]'
y = [6.57 5.15 7.15 0.18 6.48 8.95 5.24 10.54 1.24 8.05 3.56 3.4 2.18 7.16 2.32 -0.23 7.68 9.09 9.13 4.04]'
n = size(x)[1];
press = zeros(5)
X = ones(n,1)
for k=1:5
   global X = [X x.^k]
   # construct the design matrix
   P = X*((X'*X)\X')
   e = y - P*y
   press[k] = sum((e./(1 .-diag(P))).^2)
   println(press[k])
end
