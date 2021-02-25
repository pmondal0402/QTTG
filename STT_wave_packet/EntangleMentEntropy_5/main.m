clc; clear all; close all;
N = 5; % No. of elec or no. of spin.
dt = 0.1 ;
tmax = 80*dt ; 
 
% Get all possible combination of states with single flip
% 0 --> dwn
% 1 --> up
% config [0 0 0 1 0] --> [dwn dwn dwn up dwn]
tot_state = 0 ; count = 0 ;
for num_flip = 1 : N + 1 
   % elec states : Initially --> all down 
   ne = zeros(N, 1) ; ne( 1 : num_flip-1, 1 ) = 1 ;
   % spin states : Initially --> all up
   ns = ones(N, 1) ; ns( 1 : num_flip-1, 1 ) = 0 ;
 
   % elec
   config_e = 0 ; 
   config_e = get_uni_perm(ne, length(ne)) ;
   chi_f_e = get_manybodyspinor(config_e);

   % spin
   config_s = 0 ; 
   config_s = get_uni_perm(ns, length(ns)) ;
   chi_f_s = get_manybodyspinor(config_s);

   tot_state = tot_state + length(config_e(:,1))*length(config_s(:,1)) ;
   
   % Get initial manybody state
   if num_flip == 1
      psi_t0 = kron(chi_f_e, chi_f_s) ;
      psi = psi_t0 ;
      count = count + 1 ;
   else
      % get superposition of all possible states
      for ii = 1:length(config_e(:,1))
         for jj = 1:length(config_s(:,1))
            psi = psi + kron( chi_f_e(:,ii), chi_f_s(:,jj) ) ;
            count = count + 1 ;
         end
      end 
   end

end

psi = psi/sqrt(count) ;
count
% Total number of states
fprintf('Total number of states \n')
factorial(2*N)/(factorial(N)*factorial(N))
return
% Time evolution 
% fac = 1/sqrt(20) ;
% Iitialize psi
psi_t = sparse(2^(2*N),1) ;
rho = sparse(2^(2*N), 2^(2*N));
ind = 1 ; 
% Define rho
for tt = 0:dt:tmax
  tt
  % Superposition of psi starts after 10*dt time 
  if tt < 40*dt
    psi_t = psi_t0 ; 
  else
    psi_t = psi ;
  end
 rho = 0 ; 
 rho = sparse(psi_t)*sparse(psi_t') ;
 % Compute entanglement entropy between spin 1, 3
  % We need to trace out the 3 electron first
  rho_tot = rho ;
  for elec_num = 1:N
    rho_loc1 = 0 ; 
    rho_loc1 = trace_elec(rho_tot) ; 
    rho_tot = rho_loc1 ;
  end

  % Get rho for 1st local spin
  rho1_f = rho_tot ;
  for ii = 1:N-1
     rho1 = 0 ;
     rho1 = trace_elec(rho1_f) ;
     rho1_f = rho1 ;  
  end
  % Get entropy for 1st spin
  S1(ind) = get_entropy(rho1_f) ;
 
  % Trace out local spin 2 
  % rho_tot is the traced density matrix
  % Obtain rho of first and last spins
  rho_s = rho_tot ;
  for middle_s = 1:N-2
   rho_m = 0 ;
   rho_m = trace_middle_spin(rho_s) ;
   rho_s = rho_m ;
  end 
  rho_1_l =  rho_s ;%trace_middle_spin(rho_tot);
  S_1_l(ind) = get_entropy(rho_1_l) ;  
  time(ind) = tt ;

 % Get rho for last spin
 rho_l = trace_elec(rho_1_l) ;  
 Sl(ind) = get_entropy(rho_l) ;
 
 % Compute entanglement entropy for last floor(N/2) subsystem
  rho_hlf = rho_tot ;
  for half_s = 1:floor(N/2)
    rho_hlf2 = 0 ;
    rho_hlf2 = trace_last_spin(rho_hlf) ;
    rho_hlf = rho_hlf2 ;
  end 
  rho_12 = rho_hlf ; %trace_last_spin(rho_tot) ;
  S_12(ind) = get_entropy(rho_12) ;
  ind = ind + 1 ;
end

% Save data
data(:,1) = time ; data(:,2) = S1; data(:,3) = Sl; data(:,4) = S1 + Sl - S_1_l;
data(:,5) = S_1_l; data(:,6) = S_12 ;
dlmwrite(strcat('IE_Ns_', num2str(N), '.txt'), data) ;
