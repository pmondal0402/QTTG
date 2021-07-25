function [avg1 avg2] = F_fock_avg_e(rho_elec, N, Nspins)

avg1 = zeros(1,3*Nspins);
avg2 = zeros(1,Nspins);


dim_e = (4^N);

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

   pos = (a-1)*3 + 1;
   avg1(1,pos) = trace(rho_elec*Sx_e) ; 
   avg1(1,pos+1) = trace(rho_elec*Sy_e) ; 
   avg1(1,pos+2) = trace(rho_elec*Sz_e) ; 

   avg2(1,a) = trace(rho_elec*Numb) ; 

end


end

