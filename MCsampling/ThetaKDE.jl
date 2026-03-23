module ThetaKDE

export kde
using StatsBase, Random, NaNStatistics, FFTW, Optim, Plots, Distributions

function dct1d(data)
   # computes the discrete cosine transform of the column vector data
   nrows = length(data)
   # Compute weights to multiply DFT coefficients
   weight = vcat(1, 2 * (exp.(-im * (1:nrows-1) * pi / (2 * nrows))))
   # Re-order the elements of the columns of x
   data = [data[1:2:end]; data[end:-2:2]]
   # Multiply FFT by weights:
   return real(weight .* fft(data))
end

function idct1d(data)
   #computes the inverse discrete cosine transform
   nrows = length(data)
   # Compute weights
   weights = nrows * exp.(im * (0:nrows-1) * pi / (2 * nrows))
   # Compute x tilde using equation (5.93) in Jain
   data = real(ifft(weights .* data))
   #  Re-order elements of each column according to equations (5.93) and
   # (5.94) in Jain
   out = zeros(nrows)
   out[1:2:nrows] = data[1:div(nrows, 2)]
   out[2:2:nrows] = data[nrows:-1:div(nrows, 2)+1]
   return out / nrows
end

function mise(h, a2, a2k, N)
   # unbiased MISE value for a given "h"
   B = (1:length(a2)) .^ 2 * pi^2 * h^2
   val = 0.5 * sum(a2 .* (exp.(-B / 2) .- 1) .^ 2)
   val = val + 1 / (2 * N) * sum(exp.(-B) .* (2 .+ a2k .- a2))
   return val
end

function kde(data, n, MIN, MAX, h) # 5 arguments fixed bandwidth
   N = length(data)
   R = MAX - MIN
   dx = R / n
   weights, bin = histcountindices(data, MIN-0.00000001:dx:MAX+0.000000001)
   xmesh = range(MIN, MAX, length=n) .+ dx / 2

   initial_data = weights / N
   a = dct1d(initial_data) # discrete cosine transform of initial data

   a2k = a[3:2:n-1] # compute double frequency coefficient

   # compute unbiased estimator to squared coefficient
   a2 = N / (N - 1) * a[2:div(n, 2)] .^ 2 .- (2 .+ a2k) / (N - 1)

   # smooth the discrete cosine transform of initial data using "h"
   a_t = a .* exp.(-(0:n-1) .^ 2 * pi^2 * (h * R)^2 / 2)
   # now apply the inverse discrete cosine transform
   density = idct1d(a_t) / dx

   return xmesh, density, h
end

function kde(data, n, MIN, MAX) # 4 arguments automatic bandwidth 
   N = length(data)
   R = MAX - MIN
   dx = R / n
   weights, bin = histcountindices(data, MIN-0.00000001:dx:MAX+0.000000001)
   xmesh = range(MIN, MAX, length=n) .+ dx / 2

   initial_data = weights / N
   a = dct1d(initial_data) # discrete cosine transform of initial data

   a2k = a[3:2:n-1] # compute double frequency coefficient

   # compute unbiased estimator to squared coefficient
   a2 = N / (N - 1) * a[2:div(n, 2)] .^ 2 .- (2 .+ a2k) / (N - 1)

   # now compute the optimal bandwidth^2 minimizing unbiased MISE
   mise1 = h -> mise(exp(h) * R, a2, a2k, N)
   xx = -40:0.1:0
   vals = mise1.(xx)
   dv = diff(vals)
   xmin = xx[findmin(dv)[2]]
   xmax = xx[findmax(dv)[2]]
   res = optimize(mise1, xmin, xmax)
   h = exp(res.minimizer)

   # smooth the discrete cosine transform of initial data using "h"
   a_t = a .* exp.(-(0:n-1) .^ 2 * pi^2 * (h * R)^2 / 2)
   # now apply the inverse discrete cosine transform
   density = idct1d(a_t) / dx

   return xmesh, density, h
end

function kde(data, n, MIN, MAX, resamp::Bool) #automatic bandwidth, resampled data 
   if resamp
      udata = unique(data)
      uMin = minimum(udata)
      uMax = maximum(udata)
      xmesh, density, uh = kde(udata, n, uMin, uMax)
      xmesh, density, h = kde(data, n, MIN, MAX, uh)
   else
      xmesh, density, h = kde(data, n, MIN, MAX)
   end
   return xmesh, density, h
end

function kde(data, n)
   MAX = maximum(data)
   MIN = minimum(data)
   return kde(data, n, MIN, MAX)
end

function kde(data; plt=false)
   n = 2^14
   xmesh, density, h = kde(data, n)
   plt ? display(plot(xmesh, density)) : nothing
   return xmesh, density, h
end

function kde(data; plt=false, res=false)
   n = 2^14
   MIN = minimum(data)
   MAX = maximum(data)
   xmesh, density, h = kde(data, n, MIN, MAX, res)
   plt ? display(plot(xmesh, density)) : nothing
   return xmesh, density, h
end

end

# using .ThetaKDE, Plots, Statistics, Random
# Random.seed!(2)
# N = 100; K = 5000
# xorg = tan.(pi*(0.5 .- rand(N))) # original data
# medxorg = median(xorg)
# meanxorg = mean(xorg)
# x = zeros(N); mx = zeros(K)
# for i in 1:K  # resampling the data
#     ind = ceil.(Int64,N*rand(N)) # draw random indices
#     global x = xorg[ind];
#     mx[i] = median(x)
# end
# kde(x,plt=true,res=true);