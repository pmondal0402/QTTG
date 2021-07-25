function W = F_sqrt_series(A,Adag,nterms,s)
dim = length(A);
W = zeros(dim,dim);
% Add the series expansion of the operator
for j=0:nterms
    num = (factorial(2*j));
    den = ((4^j)*(factorial(j))^2*(1-2*j));
    W = W + (num/den)*mpower(Adag*A,j)/(2*s)^j;
end

end
