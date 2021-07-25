function  [blocks, V_tot, D_tot] = F_block_spin_eig(H_tot,Sz_op,N_op)
dim = length(H_tot);
% Construct a loop which finds when does the spin eigen value change
% and the number of electrons changes

V_tot = zeros(dim,dim);
D_tot = zeros(dim,dim);

spin_old = -3;
ne_old = 0;
blocks(1,1) = 1; 
b=1;
for a=1:dim
   my_state = zeros(dim,1); my_state(a,1) = 1;
   % find the spin and number of electrons in this states
   spin = my_state'*Sz_op*my_state;
   ne   = my_state'*N_op*my_state;
   
   % If the spin value changes you have reached a new block
   if ne==0 && abs(spin-spin_old)>10^(-8)
      blocks(b,2) = a-1;
      % Increment to declare next block
      b=b+1;
      % Load the value of first index of the block
      blocks(b,1) = a;
      % Load the new value of spin and number of electrons
      spin_old = spin;
      ne_old = spin;
   end   
   if ne~=0 && abs(spin-spin_old)>10^(-8)
      % Load the value of second index of the block 
      blocks(b,2) = a-1;
      % Increment to declare next block
      b=b+1;
      % Load the value of first index of the block
      blocks(b,1) = a;
      % Load the new value of spin and number of electrons
      spin_old = spin;
      ne_old = spin;
   end    
end

% Now that the blocks are made we will extract each one of them and
% diagonalize Hamiltonian in them
[x, y] = size(blocks);
tot_b = x;
blocks(tot_b,2) =dim;
for b=1:tot_b
   posi = blocks(b,1); posf = blocks(b,2);
   Hb = H_tot(posi:posf,posi:posf);
   [Vb, Db] = schur(Hb);
   V_tot(posi:posf,posi:posf) = Vb;
   D_tot(posi:posf,posi:posf) = Db;
end

% By now the eigen states must have been loaded
end