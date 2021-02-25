function [ rho_loc ] = trace_middle_spin(rho_loc)
% Traces middle spin degree of freedom
% Tracing over 2nd last spin only
dim_s = length(rho_loc) ;

for ii = 1:1 % 3
  rho_1_last = 0 ;
  rho_1_last = zeros(dim_s/(2^ii)) ;
  row_sh = 0 ;
  for i1 = 1:length(rho_1_last)/2 % number of blocks computed along columns
     col_sh = 0 ;
     for j1 = 1:length(rho_1_last)/2 % number of blocks computed along row
          rho_1_last(1 + row_sh/2:2 + row_sh/2, ...
                     1 + col_sh/2:2 + col_sh/2) =...
             rho_loc(1 + row_sh:2 + row_sh, ...
                     1 + col_sh:2 + col_sh)...
         + rho_loc(3 + row_sh:4 + row_sh,...
                   3 + col_sh:4 + col_sh) ;
         col_sh = col_sh + 4 ;
     end
     row_sh = row_sh + 4 ;
  end
   fprintf('Partial trace over 2nd last spin is done, resetting rho_loc \n')
   % reset rho_loc
   rho_loc = rho_1_last ;
end
