function [rho_loc] = first_trce(Ns, psi)
% Returns trace of all electron degree of freedom
% Ns  : No. of local spins
% psi : total manybody wavefunction

dim_s = 2^Ns ; 

% Initial psi index

rho_loc = 0 ;
for ii  = 1:dim_s
  ii
  ind_i = 0 ;
  ind_i = ii:dim_s:length(psi) ;

  for jj = 1:dim_s;
     ind_j = 0 ;
     ind_j = jj:dim_s:length(psi) ;
     rho_loc(ii, jj) = sum(psi(ind_i,1).*conj(psi(ind_j,1))) ;
  end
end

