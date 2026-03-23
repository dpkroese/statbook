function sfran_loglike(mu,sigma2_mu,sigma2,y)
   d, ni = size(y) 
   Sigmai = sigma2*diagm(ones(ni)) .+ sigma2_mu*ones(ni,ni)
   l = -(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai))
   for i=1:d
       yi = y[i,:];
       l = l - .5*(yi .- mu)'*(Sigmai\(yi .- mu))
   end
   return l
end
