function[H_ext ] =  applied_field(N_loc_spin);
tot_site = N_loc_spin;
sigmaz = [1 0;0 -1];
H_z_applied = zeros(2^tot_site,2^tot_site);
for i = 1:tot_site
last_part = kron(sigmaz,eye(2^(tot_site-i)));
%size(last_part)
H_z_applied = H_z_applied + kron(eye(2^(i-1)),last_part);
end



H_ext =H_z_applied; 
