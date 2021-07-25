function[rho_123] = get_trace_last_spin_gen(rho_1_last, dof)
% Traces over last spin dof
% dof = 3^4 for N = 4 spins with spin value S = 1 ;  
% sum diagonal elements of the matrix 

% Uncomment if rho is given instead of psi vector 
for i1  = 1:length(rho_1_last)/dof
  for j1 = 1:length(rho_1_last)/dof
    rho_123(i1, j1) = sum(diag(rho_1_last(1 + (i1-1)*dof : i1*dof, ...
                                          1 + (j1-1)*dof : j1*dof))) ; 
  end
end

