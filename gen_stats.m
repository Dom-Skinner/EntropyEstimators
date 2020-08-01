function [trans,cond_trans] = gen_stats(X,T,n)
X = mod(X,n);
X(X==0) = n;

i = 1;
while i < length(X)
    if X(i) == X(i+1)
        X(i+1) = [];
    else
        i = i+ 1;
    end
end
        
trans = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            trans(i,j) = sum((X(1:end-1) == i) & (X(2:end) == j))/T;
        end
    end
end
    
cond_trans = zeros(n,1);
cond_trans(1) = sum((X(1:end-2) == n) & (X(2:end-1) == 1) & (X(3:end) == 2))/T;
cond_trans(n) = sum((X(1:end-2) == n-1) & (X(2:end-1) == n) & (X(3:end) == 1))/T;
for i = 2:n-1
    cond_trans(i) = sum((X(1:end-2) == i-1) & (X(2:end-1) == i) & (X(3:end) == i+1))/T;
end
end
