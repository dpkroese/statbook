function y = invgampdf(z,a,l)
y = l^a*z.^(-a-1).*exp(-l./z)/gamma(a);
end