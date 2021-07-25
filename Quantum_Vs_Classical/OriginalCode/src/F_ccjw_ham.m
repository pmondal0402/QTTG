function H = F_ccjw_ham(Jsd,Nspins,N)

dim_e = (4^N);
dim_sp = 2^Nspins;
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

for a=1:N
    Sx_e(:,:,a) = 0.5*(Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sy_e(:,:,a) = 1i*0.5*(-Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sz_e(:,:,a) = 0.5*(Cdag_up(:,:,a)*C_up(:,:,a) - Cdag_dn(:,:,a)*C_dn(:,:,a));  
end

%----------------------------------------------
% Creation-Annhilation operators for HP boson
% ---------------------------------------------
% Construction of Local Spin Operators
%----------------------------------------------
% Construct the fermionic jw operators
af = zeros(dim_sp,dim_sp,Nspins);
adagf = zeros(dim_sp,dim_sp,Nspins);

for c=1:Nspins
    af(:,:,c) = (F_jw_adag(c,Nspins))';
    adagf(:,:,c) = F_jw_adag(c,Nspins);
end


% Local Spin Operators
Sx_sp = zeros(dim_sp,dim_sp,Nspins);
Sy_sp = zeros(dim_sp,dim_sp,Nspins);
Sz_sp = zeros(dim_sp,dim_sp,Nspins);

for c=1:Nspins

    Sp = F_jw_sp(a,Nspins);
    Sm = F_jw_sm(a,Nspins);
    
    Sx_sp(:,:,c) = (Sp + Sm)*0.5;
    Sy_sp(:,:,c) = (Sp - Sm)*0.5*(-1i);
    Sz_sp(:,:,c) = -(0.5*eye(dim_sp) - adagf(:,:,c)*af(:,:,c));

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
