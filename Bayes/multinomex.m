%multinomex.m
clear all
x = [53,57,147,93,38,113];
N = 10000;
p = zeros(N,2,3); a = zeros(N,2,3);
alpha = x + 1;
for i=1:N
    r = dirichrnd(alpha);
    h = [r,1-sum(r)];
    p(i,:,:) = reshape(h',3,2)';
end
for i=1:2
    p_row(i,:) = sum(p(:,i,:),3); 
end
for j=1:3
    p_col(j,:) = sum(p(:,:,j),2);
end
for i=1:2
    for j=1:3
        a(:,i,j) = p(:,i,j) - p_row(i,:)'.*p_col(j,:)';
    end
end
hold on
for j=1:3
    kde(a(:,1,j));
end