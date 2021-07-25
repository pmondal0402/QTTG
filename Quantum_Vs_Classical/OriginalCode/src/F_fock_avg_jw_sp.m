function [avg magno] = F_fock_avg_2by2_sp(psi_t,N,Nspins)

avg = zeros(1,3*Nspins);
dim_e = (4^N);
dim_sp = 2^Nspins;
dim = dim_e*dim_sp;

%----------------------------------------------
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


magno = 0;
for c=1:Nspins

    Sp = F_jw_sp(c,Nspins);
    Sm = F_jw_sm(c,Nspins);

    Sx_sp = (Sp + Sm)/2;
    Sy_sp = (Sp - Sm)/(2*1i);
    Sz_sp = -(0.5*eye(dim_sp) - adagf(:,:,c)*af(:,:,c));

    Sx_sp = kron(eye(dim_e),Sx_sp);
    Sy_sp = kron(eye(dim_e),Sy_sp);
    Sz_sp = kron(eye(dim_e),Sz_sp);

   pos = (c-1)*3 + 1;
   avg(1,pos) = psi_t'*Sx_sp*psi_t;
   avg(1,pos+1) = psi_t'*Sy_sp*psi_t;
   avg(1,pos+2) = psi_t'*Sz_sp*psi_t;

    
end


end
