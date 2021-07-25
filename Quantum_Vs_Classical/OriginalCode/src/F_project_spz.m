function Pspz = F_project_spz(Nspins,N,tot_sz,n,s)
% Construct the total Sz operator
Sz_op  = F_tot_sz_op(Nspins,N,n,s);
dim = length(Sz_op);
% find the eigen states of this operators
[V D] = eig(Sz_op);
% Now construct the Projection operator
Pspz = zeros(dim,dim);

for a=1:dim
    if D(a,a)==tot_sz
       vec = V(:,a);
       Pspz = Pspz + vec*vec'; 
    end
end

end

