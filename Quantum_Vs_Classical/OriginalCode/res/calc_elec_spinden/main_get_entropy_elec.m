% -----------------------------------------------------------------------                                                                                                          
% Note : Here I am calculating the following :
% You can plot both entropy of whole electronic density matrix, its spin 
% density matrix, and |P| from its spin density matrix. 
% -----------------------------------------------------------------------
% TODO : entropy of whole electronic density matrix
% |P| from its spin density matrix 

clc; clear all; close all;
% Number of electrons
N = 4 ; 
e_dof = 4 ; 
s_dof = 3 ; 
file_sp = fopen('../elec_spin.txt','w') ; 
file_sp2 = fopen('../elec_pop.txt','w') ; 
Nspins = N ;
s = 1 ;  
n = 2*s + 1 ; 
 
% Read data
for sv_frm = 1:20001% 20001 % TODO change this val to 20001 once job is complete
  sv_frm
  clear S ;
  path = sprintf('../saved_rho/psi_t%03d.mat', sv_frm) ; 
  % ../saved_rho/rho_t%03d.mat', sv_frm) ;
  S = load(path) ;
  psi_t = S.psi ; 

  rho_s = psi_t*sparse(psi_t');
  % Compute elec spin density 
  % Get elec density matrix 
  rho_elec = get_trace_last_spin_gen(rho_s, s_dof^N) ;

  [avg1 avg2] = F_fock_avg_e(rho_elec, N, Nspins) ; 

  % Save data 
  fprintf(file_sp,...
'%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n', ...
avg1);   
  fprintf(file_sp2, '%22.16f %22.16f %22.16f %22.16f\n', avg2) ; 
end

% NOTE : Does the electron reach 4thr site at t = 40 fs and reflects back around
% 80 fs and that's how we see pattern in entanglement spin sub ? 
'Job is done ! '
