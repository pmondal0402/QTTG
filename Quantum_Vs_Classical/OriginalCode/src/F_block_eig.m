function [V_tot, D_tot, blocks] = F_num_block_eig(H_tot,Dn,N)
% The dimension of the Hamiltoninan
dim = length(H_tot);
% The output eigen states and eigen values
V_tot = zeros(dim,dim);
D_tot = zeros(dim,dim);

% Find the blocks one by one
blocks = zeros(2*N+1,2);

for number=0:2*N
   % Take Dn and scan to see when its value deviates from the current particle number
   count=0;
   count2 = 0;
   for scan=1:dim
       if number==0
          % Find when the scan changes value
          blocks(number+1,1)=1;
          if Dn(scan,scan)~=0 && count==0
              blocks(number+1,2)= scan-1;
              count=count+1;
          end
       end
       
       
       if number>0 && number<2*N
           % Find where does the scan value change
           blocks(number+1,1) = blocks(number,2)+1;
           if Dn(scan,scan)~=number &&...
              scan>blocks(number+1,1) && ...
              count2==0
              blocks(number+1,2)= scan-1;
              count2=count+1;
           end
       end
       
       if number==2*N
          % Find where does the scan value change
           blocks(number+1,1) = blocks(number,2)+1;
           blocks(number+1,2) = dim;
       end
   end
   
   % Take our Hamiltonian and block diagonalize it
   posi = blocks(number+1,1); posf = blocks(number+1,2);
   Hb = H_tot(posi:posf,posi:posf);
   [Vb,Db] = eig(Hb);
   V_tot(posi:posf,posi:posf) = Vb;
   D_tot(posi:posf,posi:posf) = Db;
end


end

