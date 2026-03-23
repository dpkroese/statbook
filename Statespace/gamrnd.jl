function gamrnd(a,c) #c is scale not rate
  n = length(c)
  x = zeros(n)
  for i=1:n
    x[i] = rand(Gamma(a,c[i]))
  end
  return x
end
