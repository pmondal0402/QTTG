function [Hamil_heis] = H_heisenbergp(Lx,Ly,N_loc_spin,Jx,Jy,Jz, pbc)


% Name the file same as function i.. H_sd.m in this case
% Lx,Ly,N_loc_spin,pos_loc_i a input variables
% Lx,Ly = latice size
%N_loc_spin = no. of local spins involved
% pos_loc_i = position of first local spin 
% Hamil_sd = Output variable

  sigma_x = [0 1;1 0];
  sigma_y = [0 -1i;1i 0];
  sigma_z = [1 0;0 -1];
  spin_conf = 2^(N_loc_spin);

%c = Lx*Ly;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%          Calculating spin Hamiltonian first

% .................................................................................................................

%ixx = [1:N_loc_spin];

% Calculating horizontal interactions first ....................

H_heisenberg = sparse(spin_conf,spin_conf);

% Calculating interaction terms S1.S2, S2.S3,S3.S4......S5.S6,S6,S7,S7.S8,S9.S10

x_12 = kron(sparse(sigma_x),sigma_x);
y_12 = kron(sparse(sigma_y),sigma_y);
z_12 = kron(sparse(sigma_z),sigma_z);

	
for i = 1:(N_loc_spin -1 )

	if i == Lx	% Ignoring term S5.S6 as they are not nearst neighbours.............use mod(m,n) function for general case  
		continue
	end

	last_part = ((N_loc_spin - 2) - (i-1));

	x = kron(x_12,eye(2^last_part));
	y = kron(y_12,eye(2^last_part));
	z = kron(z_12,eye(2^last_part));

    
	H_heisenberg = H_heisenberg + Jx*kron(eye(2^(i-1)),x) + Jy*kron(eye(2^(i-1)),y) + Jz*kron(eye(2^(i-1)),z);

end



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if Ly ~= 1

% Calculating vertical interactions  ....................
% Calculating terms S1.S6, S2.S6, S3.S7,S4.S8
Jz

x_15 = kron(sparse(sigma_x),eye(2^(Lx-1)) );
y_15 = kron(sparse(sigma_y),eye(2^(Lx-1)) );
z_15 = kron(sparse(sigma_z),eye(2^(Lx-1)) );



x_15 = kron(x_15,sigma_x);
y_15 = kron(y_15,sigma_y);
z_15 = kron(z_15,sigma_z);


	for i  = 1: (N_loc_spin/Ly)

	 last_partt = ((Lx-1) - (i-1));	
	 x1 = kron(x_15,eye(2^last_partt));
	 y1 = kron(y_15,eye(2^last_partt));
	 z1 = kron(z_15,eye(2^last_partt));

		%a1 = x1 +y1 +z1;

		H_heisenberg = H_heisenberg + Jx*kron(eye(2^(i-1)),x1) + Jy*kron(eye(2^(i-1)),y1) + Jz*kron(eye(2^(i-1)),z1);

	end

end
pbc

if pbc >= 1
  % Add periodic boundary condition
  S1x = kron(eye(2^(N_loc_spin - 1)), sigma_x) ;  
  S1y = kron(eye(2^(N_loc_spin - 1)), sigma_y) ;
  S1z = kron(eye(2^(N_loc_spin - 1)), sigma_z) ;
  
  SNx = kron(sigma_x, eye(2^(N_loc_spin - 1))) ; 
  SNy = kron(sigma_y, eye(2^(N_loc_spin - 1))) ; 
  SNz = kron(sigma_z, eye(2^(N_loc_spin - 1))) ; 
  
  Hamil_heis = H_heisenberg + Jx*S1x*SNx ...
                            + Jy*S1y*SNy ...
                            + Jz*S1z*SNz ;
else
  Hamil_heis = H_heisenberg ; 
end

% Hamil_heis = H_heisenberg ; 

