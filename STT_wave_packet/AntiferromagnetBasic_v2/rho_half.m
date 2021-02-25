function[rho_hlf] = rho_half(rho_tot, N)

% Traces out half the degree of freedom of input rho_tot
rho_hlf = rho_tot ;
for half_s = 1 : floor(N/2)
  rho_hlf2 = 0 ; 
  rho_hlf2 = trace_last_spin(rho_hlf) ;
  rho_hlf = rho_hlf2 ; 
end

