function H = F_ccdag_ham2(gamma,N)
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
H  = zeros(dim,dim);
H1 = zeros(N,N); 

% Single particle noninteracting Hamiltonian
for a=1:N-1
 H1(a,a+1) = -gamma;
 H1(a+1,a) = -gamma;
end

D = eig(H1);

% Set up the non-interacting Hamiltonian in energy basis
for a=1:N
     H = H + D(a,1)*(Cdag_up(:,:,a)*C_up(:,:,a))  ...
           + D(a,1)*(Cdag_dn(:,:,a)*C_dn(:,:,a));
end

end
