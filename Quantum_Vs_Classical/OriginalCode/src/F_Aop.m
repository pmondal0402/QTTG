function A =F_Aop(n)
A = zeros(n,n);
for c=1:n-1
    A(c,c+1) = sqrt(c);
end
end
