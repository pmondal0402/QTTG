function [avg magno] = F_fock_avg_2by2_sp(psi_t,N,Nspins)

avg = zeros(1,3*Nspins);
dim_e = (4^N);
dim_sp = 2^Nspins;
dim = dim_e*dim_sp;

%----------------------------------------------
% Creation-Annhilation operators for HP boson
% ---------------------------------------------
% Construction of Local Spin Operators
%----------------------------------------------

magno = 0;
for c=1:Nspins

    Sx_sp = F_2by2_SxSySz(c,Nspins,1);
    Sy_sp = F_2by2_SxSySz(c,Nspins,2);
    Sz_sp = F_2by2_SxSySz(c,Nspins,3);

    Sx_sp = kron(eye(dim_e),Sx_sp);
    Sy_sp = kron(eye(dim_e),Sy_sp);
    Sz_sp = kron(eye(dim_e),Sz_sp);
    mag = kron(eye(dim_e),Sz_sp);

   pos = (c-1)*3 + 1;
   avg(1,pos) = psi_t'*Sx_sp*psi_t;
   avg(1,pos+1) = psi_t'*Sy_sp*psi_t;
   avg(1,pos+2) = psi_t'*Sz_sp*psi_t;

    
end


end
