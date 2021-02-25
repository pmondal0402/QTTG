function [Hamil_heis] = H_heisenberg(Lx,Ly,N_loc_spin,Jx,Jy,Jz)


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



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% ..............................................................

% Calculate vertical interaction terms 


% .............................................................

Hamil_heis = H_heisenberg;
%Hamil_heis =kron(eye(Lx*Ly), H_heisenberg);%H_heisenberg;


%Hamil_heis = H_heisenberg;




% Ignore rest of code for now 
% .................................................................................................................


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%{






  pos_loc_f = pos_loc_i+N_loc_spin-1;      % Position where local spin stops
  
% Calculating local spins in total spin space 

  for i = 1:N_loc_spin
        [spin_loc_x{i,1:(2^(i+1))}] = deal(kron(eye((2^(i))),sparse(sigma_x))) ;
        [spin_loc_y{i,1:(2^(i+1))}] = deal(kron(eye((2^(i))),sparse(sigma_y))) ;
        [spin_loc_z{i,1:(2^(i+1))}] = deal(kron(eye((2^(i))),sparse(sigma_z))) ;
 
 
        % it assigns spin_loc(N_loc_spin + 1,2,3...) = default 24x24
 
             for pos = pos_loc_i:(pos_loc_f-i)
               %  spin_loc3_x = kron(spin_loc3_x,eye(2));
                % pos
                 spin_loc_x{i} = kron(spin_loc_x{i},sparse(eye(2)));
                 spin_loc_y{i} = kron(spin_loc_y{i},sparse(eye(2)));
                 spin_loc_z{i} = kron(spin_loc_z{i},sparse(eye(2)));
 
             end
 
   end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
                    Hheis = sparse(spin_conf,spin_conf);
            %...........................................

        % Considering only nearest neighbour interactions for local spin coupling 

        % Calculating nearest neighbour along vertical and horizontal direction 
	

	if Lx ~=1
		for 
                for i = 1: (N_loc_spin/Ly)-1    
                    Hheis = Hheis + spin_loc_x{i}*spin_loc_x{i+1} ...
                                        + spin_loc_y{i}*spin_loc_y{i+1} ...
                                        + spin_loc_z{i}*spin_loc_z{i+1};
                end
	end      
%............adding the boundary terms for a circular ring whih also corresponds to infinite long chain of local spins .............

		if circular == 1 && Ly == 1

			Hheis = Hheis +  spin_loc_x{N_loc_spin}*spin_loc_x{1} ...
				      +  spin_loc_y{N_loc_spin}*spin_loc_y{1} ...
				      +  spin_loc_z{N_loc_spin}*spin_loc_z{1}; 
		circular
		end

%...................................................................................................................................


	  if Ly ~= 1
              if N_loc_spin > Ly

                for i = 1:N_loc_spin -  Ly
                    Hheis = Hheis + spin_loc_x{i}*spin_loc_x{i + Ly} ...
                                        + spin_loc_y{i}*spin_loc_y{i + Ly} ...
                                        + spin_loc_z{i}*spin_loc_z{i + Ly} ;
                end

              end
         end

		if circular == 1
			Hamil_heis = (kron(eye(Lx*Ly),Hheis) + ctranspose(kron(eye(Lx*Ly),Hheis)))/4;
		else
                	Hamil_heis = (kron(eye(Lx*Ly),Hheis));% + ctranspose(kron(eye(Lx*Ly),Hheis)))/2;
		end
%}


%Hamil_heis = ((Hheis) + ctranspose((Hheis)))/2;
