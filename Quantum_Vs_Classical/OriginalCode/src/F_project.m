function P = F_project(num_p,N)
% Dimension of the Hilbert space
dim = 4^N;
% Construct the number operator
nop = F_num_op(N);
% Find the eigen states with eigen value num_p3 and take the projector
[V,D] = eig(nop);

% Initalize projection operator
P = zeros(dim,dim);
for a=1:dim
    test = real(D(a,a));
    if test==num_p
         vec = V(:,a);
         P = P + vec*vec';
    end    
end

end
