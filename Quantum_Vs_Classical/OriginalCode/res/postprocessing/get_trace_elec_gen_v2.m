function [rho_loc] = get_trace_elec_gen_v2(psi, dof)
% dof : degree of freedom with val (2*s + 1) 
% Returns trace of 1st (electron) degree of freedom
dim_s = length(psi)/dof ;

rho_loc = sparse(dim_s, dim_s) ; 
for i = 1:dof
  rho_loc = rho_loc ...
          + psi(1 + (i-1)*dim_s: i*dim_s, 1)*psi(1 + (i-1)*dim_s: i*dim_s, 1)' ; 
          % + rho(1 + (i-1)*dim_s: i*dim_s, 1 + (i-1)*dim_s: i*dim_s) ; 
end
 
