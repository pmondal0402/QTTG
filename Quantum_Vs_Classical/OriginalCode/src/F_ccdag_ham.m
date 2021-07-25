function H = F_ccdag_ham(gamma,N)
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

% Create the Hamiltonian
H = zeros(dim,dim);
for a=1:N-1
    H = H - gamma*(Cdag_up(:,:,a+1)*C_up(:,:,a) + Cdag_up(:,:,a)*C_up(:,:,a+1)) ...
          - gamma*(Cdag_dn(:,:,a+1)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_dn(:,:,a+1));
end


end
