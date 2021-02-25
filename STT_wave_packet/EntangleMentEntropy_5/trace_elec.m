function [rho_loc] = trace_elec(rho)
% Returns trace of 1st (electron) degree of freedom
dim_s = length(rho)/2 ;  
rho_loc = rho(1:dim_s, 1:dim_s) ...
        + rho(1 + dim_s: 2*dim_s, 1 + dim_s: 2*dim_s);

