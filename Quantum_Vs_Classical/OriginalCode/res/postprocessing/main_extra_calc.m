% -----------------------------------------------------------------------
% Note : Here I am calculating the following :
% You can plot both entropy of whole electronic density matrix, its spin 
% density matrix, and |P| from its spin density matrix. 
% -----------------------------------------------------------------------
% TODO : entropy of whole spin density matrix

clc; clear all; close all;
% Number of electrons
N = 4 ; 
e_dof = 4 ; 
s_dof = 3 ; 
dim1 = s_dof^(N/2) ;
dim2 = dim1 ; 
file_sp = fopen('../entropy_spin.txt','w') ;
 
% Read data
for sv_frm = 1:20001% 20001 % TODO change this val to 20001 once job is complete
  sv_frm
  clear S ;
  path = sprintf('../saved_rho/psi_t%03d.mat', sv_frm) ; 
  % ../saved_rho/rho_t%03d.mat', sv_frm) ;
  S = load(path) ;
  psi_t = S.psi ; 

  % Get rho for spins by tracing over electrons
  rho_s = get_trace_elec_gen_v2(psi_t, e_dof) ; 
  for sp = 1:N-1
    rho_m = 0 ;
    rho_m = get_trace_elec_gen(rho_s, e_dof) ; 
    rho_s = rho_m ; 
  end
  
  % save rho spin 
  % rho_sp = rho_s ; 

  % Get entropy for the same
  S_en(sv_frm) = get_entropy(rho_s) ;  
  % Save data
  data = [real(sv_frm); real(S_en(sv_frm))] 
  fprintf(file_sp,'%22.16f %22.16f\n', data); 
end

