function [avg magno] = F_fock_avg_sp(psi_t,N,Nspins,n,s,nterms)

avg = zeros(1,3*Nspins);
dim_e = (4^N);
dim_sp = n^Nspins;
dim = dim_e*dim_sp;

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

magno = 0;
for c=1:Nspins

    Spi = sqrt(2*s)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s)*af(:,:,c);
    Smi = sqrt(2*s)*adagf(:,:,c)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s);
    mag = adagf(:,:,c)*af(:,:,c);

    Sx_sp = 0.5*(Spi+Smi);
    Sy_sp = -1i*0.5*(Spi-Smi);
    Sz_sp = (s*eye(dim_sp) - adagf(:,:,c)*af(:,:,c));

    Sx_sp = kron(eye(dim_e),Sx_sp);
    Sy_sp = kron(eye(dim_e),Sy_sp);
    Sz_sp = kron(eye(dim_e),Sz_sp);
    mag = kron(eye(dim_e),mag);

   pos = (c-1)*3 + 1;
   avg(1,pos) = psi_t'*Sx_sp*psi_t;
   avg(1,pos+1) = psi_t'*Sy_sp*psi_t;
   avg(1,pos+2) = psi_t'*Sz_sp*psi_t;

   magno = magno + psi_t'*mag*psi_t;
    
end


end
