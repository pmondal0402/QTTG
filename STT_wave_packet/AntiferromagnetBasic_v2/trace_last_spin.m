function[rho_123] = trace_last_spin(rho_1_last)
% Traces over last spin dof
row_sh = 0 ;
for i1  = 1:length(rho_1_last)/2
  col_sh = 0 ;
  for j1 = 1:length(rho_1_last)/2
    rho_123(i1, j1) = rho_1_last(i1 + row_sh/2, j1 + col_sh/2) ...
                    + rho_1_last(i1 + row_sh/2 + 1, j1 + col_sh/2 + 1) ;
    col_sh = col_sh + 2 ;
  end
  row_sh = row_sh  + 2 ;
end

