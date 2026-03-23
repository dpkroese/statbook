include("../../Julia/KDE/ThetaKDE/ThetaKDE.jl")
using Distributions, .ThetaKDE, Plots
x = [53,57,147,93,38,113];
N = 10000;
p = zeros(N,2,3); a = zeros(N,2,3);
p_row = zeros(2,N); p_col=zeros(3,N);
alpha = x .+ 1;
for i in 1:N
    h = rand(Dirichlet(alpha));
    p[i,:,:] = reshape(h',3,2)';
end
for i in 1:2
    p_row[i,:] = sum(p[:,i,:],dims=2); 
end
for j in 1:3
    p_col[j,:] = sum(p[:,:,j],dims=2);
end
for k in 1:N
    for i in 1:2
        for j in 1:3
            a[k,i,j] = p[k,i,j] - p_row[i,k]*p_col[j,k];
        end
    end
end
p = plot()
for j in 1:3
    xmesh, density, h = kde(a[:,1,j])
    p = plot!(xmesh,density)
    display(p)
end