function [avg1 avg2] = F_fock_avg_e(psi_t,N,Nspins,n,s);

avg1 = zeros(1,3*Nspins);
avg2 = zeros(1,Nspins);


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

for a=1:N
    Sx_e = 0.5*(Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sy_e = 1i*0.5*(-Cdag_up(:,:,a)*C_dn(:,:,a) + Cdag_dn(:,:,a)*C_up(:,:,a));
    Sz_e = 0.5*(Cdag_up(:,:,a)*C_up(:,:,a) - Cdag_dn(:,:,a)*C_dn(:,:,a));
    Numb = Cdag_up(:,:,a)*C_up(:,:,a) + Cdag_dn(:,:,a)*C_dn(:,:,a);

    Sx_e = kron(Sx_e,eye(dim_sp));
    Sy_e = kron(Sy_e,eye(dim_sp));
    Sz_e = kron(Sz_e,eye(dim_sp));
    Numb = kron(Numb,eye(dim_sp));

   pos = (a-1)*3 + 1;
   avg1(1,pos) = psi_t'*Sx_e*psi_t;
   avg1(1,pos+1) = psi_t'*Sy_e*psi_t;
   avg1(1,pos+2) = psi_t'*Sz_e*psi_t;

   avg2(1,a) = psi_t'*Numb*psi_t;

end


end

