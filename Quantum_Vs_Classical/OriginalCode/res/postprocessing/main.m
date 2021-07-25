clc; clear all; close all;
% Number of electrons
N = 4 ; 
e_dof = 4 ; 
s_dof = 3 ; 
dim1 = s_dof^(N/2) ;
dim2 = dim1 ; 
 
% Read data
for sv_frm = 1:20001 % TODO change this val to 20001 once job is complete
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
  rho_sp = rho_s ; 

  % Get rho for half of spin system 
  for sp = 1:N-2
    rho_m2 = 0 ;
    rho_m2 = get_trace_elec_gen(rho_s, s_dof) ; 
    rho_s = rho_m2 ; 
  end
  
  % Get entropy for the same
  S_en(sv_frm) = get_entropy(rho_s) ;  
  
  % Compute Negativity of subsystem
  [neg(sv_frm), log_neg(sv_frm)] = ...
                   negativity(full(rho_sp), [dim1, dim2]) ; 
   
end

dlmwrite('../S_enHalf.txt', real(S_en)) ;
data_neg(:,1) = neg ; data_neg(:, 2) = log_neg ; 
dlmwrite('../negativity.txt', data_neg) ;  
