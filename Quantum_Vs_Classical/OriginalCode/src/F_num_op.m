function nop = F_num_op(N)
dim = 4^N;

% Creation-Annhilation operators for electrons
Cdag_up = zeros(dim,dim,N);
Cdag_dn = zeros(dim,dim,N);
C_up = zeros(dim,dim,N);
C_dn = zeros(dim,dim,N);

% Fill the operators
for a=1:N
   Cdag_up(:,:,a) = F_Cdag_up(a,N);
   Cdag_dn(:,:,a) = F_Cdag_dn(a,N);
   C_up(:,:,a) = F_Cdag_up(a,N)';
   C_dn(:,:,a) = F_Cdag_dn(a,N)';
end

% Create the number operator
nop = zeros(dim,dim);
for a=1:N
    nop = nop + (Cdag_up(:,:,a)*C_up(:,:,a) + Cdag_dn(:,:,a)*C_dn(:,:,a));
end

end
