function Stot = F_tot_sz_op(Nspins,N,n,s)

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

for a=1:N
    Sx_e(:,:,a) = 0.5*(Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sy_e(:,:,a) = 1i*0.5*(-Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sz_e(:,:,a) = 0.5*(Cdag_up(:,:,a)*C_up(:,:,a) - Cdag_dn(:,:,a)*C_dn(:,:,a));
end

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

Sz_sp = zeros(dim_sp,dim_sp,Nspins);
for c=1:Nspins
    Sz_sp(:,:,c) = (s*eye(dim_sp) - adagf(:,:,c)*af(:,:,c));
end

%------------------------
% Create the total operator
%------------------------
Stot = zeros(dim,dim);

for a=1:N
    Stot = Stot + kron(Sz_e(:,:,a),eye(dim_sp)) + kron(eye(dim_e),Sz_sp(:,:,a));
end

