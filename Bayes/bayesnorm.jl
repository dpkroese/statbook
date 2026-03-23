include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Distributions, .ThetaKDE, Plots
n = 10; 
X = randn(n); #generate the data
sample_mean = mean(X);
sample_var = var(X);
sig2 = var(X); mu = sample_mean; #initial state
N = 10^5; #sample size for Gibbs sampler 
gibbs_sample = zeros(N,2);
for k in 1:N    
    global mu = sample_mean + sqrt(sig2/n)*randn();  #draw mu
    V = sum((X .- mu).^2)/2;
    global sig2 = 1/rand(Gamma(n/2,1/V)); #draw sigma^2
    gibbs_sample[k,:] = [mu sig2];  
end
p1 = xmesh,density,bw = kde(gibbs_sample[:,1]); 
p2 = xmesh,density,bw = kde(gibbs_sample[:,2]);
plot([p1,p2],layout=(1,2))
