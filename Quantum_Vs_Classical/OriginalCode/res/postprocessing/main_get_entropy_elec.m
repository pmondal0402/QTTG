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
file_sp = fopen('../entropy_elec.txt','w') ; 
 
% Read data
for sv_frm = 1:20001% 20001 % TODO change this val to 20001 once job is complete
  sv_frm
  clear S ;
  path = sprintf('../saved_rho/psi_t%03d.mat', sv_frm) ; 
  % ../saved_rho/rho_t%03d.mat', sv_frm) ;
  S = load(path) ;
  psi_t = S.psi ; 

  rho_s = psi_t*sparse(psi_t'); 
  
  % Get rho of electron spin system 
  rho_m2 = get_trace_last_spin_gen(rho_s, s_dof^N) ; 

  % Get entropy for the same
  S_en(sv_frm) = get_entropy(rho_m2) ;  
  % Save data 
  data = [real(sv_frm); real(S_en(sv_frm))]
  fprintf(file_sp,'%22.16f %22.16f\n', data);   
end

