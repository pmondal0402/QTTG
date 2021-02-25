function [rho_mid] = rho_loc(rho_gs, ns, Ntot)
% rho_gs --> total spin density matrix
% ns --> What is the spin number in the chain
% Ntot --> Total number of spins in the chain
% Example described to compute rho2 in a spin chain of 1--2--3--4
% trace over 1
rho_mid = rho_gs ;  
      
num_frst = ns - 1 ;  % no. of first dof to be traced
for frst = 1:num_frst
    rho_mid1 = 0 ; 
    rho_mid1 = trace_elec(rho_mid) ; 
    rho_mid = rho_mid1 ; 
end
                                 
% trace over 3 and 4
for lst = 1:Ntot-ns
    rho_mid2 = 0 ; 
    rho_mid2 = trace_last_spin(rho_mid) ;
    rho_mid = rho_mid2 ; 
end
% ns
% rho_mid
