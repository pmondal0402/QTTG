function [chi_f] = get_manybodyspinor(config)
% Returns manybody spinor state
% config [0 0 0 1 0] --> [dwn dwn dwn up dwn]
% Initialize
% No. of spins
num_s = length(config(1,:)) ;
% No. of configurations
num_config = length(config(:,1)) ;
 
chi_f = zeros(2^num_s, num_config) ;
% Do for each config
for in = 1:num_config
  chi0 = 1 ;
  for jn = 1:num_s
      if config(in, jn) == 1
        chi = [1 ; 0] ; 
      else
        chi = [0; 1] ;
      end
      chi0 = kron(chi0, chi) ;
  end
  % Now chi0 is the manybody spinor state
  %size(chi0)
  %size(chi_f)
  chi_f(:,in) = chi0 ;
end

