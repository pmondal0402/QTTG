function H = F_cacadag_ham2(Jsd,Nspins,N,n,s)

dim_e = (4^N);
dim_sp = n^Nspins;
dim = dim_e*dim_sp;

%----------------------------------------------
% Creation-Annhilation operators for electrons
% ---------------------------------------------
% Construction of Electronic Spin Operators
%----------------------------------------------
Cdag_up = zeros(dim_e,dim_e,N);
Cdag_dn = zeros(dim_e,dim_e,N);
C_up = zeros(dim_e,dim_e,N);
C_dn = zeros(dim_e,dim_e,N);

% Fill the operators
for a=1:N
   Cdag_up(:,:,a) = F_Cdag_up(a,N);
   Cdag_dn(:,:,a) = F_Cdag_dn(a,N);
   C_up(:,:,a) = F_Cdag_up(a,N)';
   C_dn(:,:,a) = F_Cdag_dn(a,N)';
end

% Electronic Spin operators
Sx_e = zeros(dim_e,dim_e,N);
Sy_e = zeros(dim_e,dim_e,N);
Sz_e = zeros(dim_e,dim_e,N);

% Create the electronic spin operators
for a=1:N
    Sx_e(:,:,a) = F_e_spin_op(a,C_up,C_dn,N,1);
    Sy_e(:,:,a) = F_e_spin_op(a,C_up,C_dn,N,2);
    Sz_e(:,:,a) = F_e_spin_op(a,C_up,C_dn,N,3);  
end

%['Writing files ...']
%dlmwrite('../testing/sx_e_site1_N3_utk.txt',kron(Sx_e(:,:,1),eye(dim_sp)));
%dlmwrite('../testing/sy_e_site1_N3_utk.txt',kron(Sy_e(:,:,1),eye(dim_sp)));
%dlmwrite('../testing/sz_e_site1_N3_utk.txt',kron(Sz_e(:,:,1),eye(dim_sp)));


%----------------------------------------------
% Creation-Annhilation operators for HP boson
% ---------------------------------------------
% Construction of Local Spin Operators
%----------------------------------------------

% Construct the ladder operators
af = zeros(dim_sp,dim_sp,Nspins);
adagf = zeros(dim_sp,dim_sp,Nspins);

% Fill the operators
for a=1:N
    af(:,:,a) = F_a_fock(a,Nspins,n);
    adagf(:,:,a) = (F_a_fock(a,Nspins,n))';
end

% Local Spin Operators
Sx_sp = zeros(dim_sp,dim_sp,Nspins);
Sy_sp = zeros(dim_sp,dim_sp,Nspins);
Sz_sp = zeros(dim_sp,dim_sp,Nspins);

for c=1:Nspins

    Spi = sqrt(2*s)*sqrtm(eye(dim_sp) - adagf(:,:,c)*af(:,:,c)/(2*s))*af(:,:,c);
    Smi = sqrt(2*s)*adagf(:,:,c)*sqrtm(eye(dim_sp) - adagf(:,:,c)*af(:,:,c)/(2*s));

    Sx_sp(:,:,c) = 0.5*(Spi+Smi);
    Sy_sp(:,:,c) = -1i*0.5*(Spi-Smi);
    Sz_sp(:,:,c) = (s*eye(dim_sp) - adagf(:,:,c)*af(:,:,c));

end

%------------------------
% Create the Hamiltonian
%------------------------
H = zeros(dim,dim);

for a=1:N
    H = H - Jsd*kron(Sx_e(:,:,a),Sx_sp(:,:,a)) ...
          - Jsd*kron(Sy_e(:,:,a),Sy_sp(:,:,a)) ...
          - Jsd*kron(Sz_e(:,:,a),Sz_sp(:,:,a));
end

end
