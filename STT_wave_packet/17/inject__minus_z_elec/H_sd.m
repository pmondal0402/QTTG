function [Hamil_sd] = H_sd(Lx,Ly,N_loc_spin,pos_loc_i)
% Name the file same as function i.. H_sd.m in this case
% Lx,Ly,N_loc_spin,pos_loc_i a input variables
% Lx,Ly = latice size
%N_loc_spin = no. of local spins involved
% pos_loc_i = position of first local spin 
% Hamil_sd = Output variable

  sigma_x = [0 1;1 0];
  sigma_y = [0 -1i;1i 0];
  sigma_z = [1 0;0 -1];
  spin_conf = 2^(N_loc_spin + 1);

  pos_loc_f = pos_loc_i+N_loc_spin-1;      % Position where local spin stops
  
% Calculating local spins in total spin space 

  for i = 1:N_loc_spin
        [spin_loc_x{i,1:(2^(i+1))}] = deal(kron(eye((2^i)),sparse(sigma_x))) ;
        [spin_loc_y{i,1:(2^(i+1))}] = deal(kron(eye((2^i)),sparse(sigma_y))) ;
        [spin_loc_z{i,1:(2^(i+1))}] = deal(kron(eye((2^i)),sparse(sigma_z))) ;
 
 
        % it assigns spin_loc(N_loc_spin + 1,2,3...) = default 24x24
 
             for pos = pos_loc_i:(pos_loc_f-i)
               %  spin_loc3_x = kron(spin_loc3_x,eye(2));
                % pos
                 spin_loc_x{i} = kron(spin_loc_x{i},sparse(eye(2)));
                 spin_loc_y{i} = kron(spin_loc_y{i},sparse(eye(2)));
                 spin_loc_z{i} = kron(spin_loc_z{i},sparse(eye(2)));
 
             end
 
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Calculating Electron spins in total spin space
 
                 sigma_e_x = kron(sparse(sigma_x),eye(2^N_loc_spin));
                 sigma_e_y = kron(sparse(sigma_y),eye(2^N_loc_spin));
                 sigma_e_z = kron(sparse(sigma_z),eye(2^N_loc_spin));

 %...................Hsd coupling ...........................................
                 site_pos = sparse(Lx*Ly,1);
                 Hsd =  sparse( Lx*Ly*spin_conf,Lx*Ly*spin_conf );
 
                 for i = 1:N_loc_spin
                 site_pos(pos_loc_i +i-1,1) = 1;
                 %full(site)
                 Hsd = Hsd + kron((site_pos*site_pos'),(sigma_e_x*spin_loc_x{i} + sigma_e_y*spin_loc_y{i} + sigma_e_z*spin_loc_z{i}));
                 site_pos(pos_loc_i +i-1,1) = 0;
                 end
                 Hamil_sd = (Hsd); %+ ctranspose(Hsd))/2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

