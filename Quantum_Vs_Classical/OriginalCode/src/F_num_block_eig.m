function [blocks, U2] = F_num_block_eig(H_tot,Dn,N,Sz_op)
% The dimension of the Hamiltoninan
dim = length(H_tot);

U2 = zeros(dim,dim);
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
   Sb = Sz_op(posi:posf,posi:posf);
   
   % In order to Also make them eigen states of spin operator
   % we must do unitary transformation on Hb using block U2
   
   [Vs, Ds] = eig(Sb); % Block of U2 is Vs  
  
   % Unitary Transformation of second type
   U2(posi:posf,posi:posf) = Vs;
   
%    % Hamiltonian in the block
%    Hb = H_tot(posi:posf,posi:posf);
%    
%    % Hamiltonian further blocked by spin
%    Hbb = Vs'*Hb*Vs;
%    Sbb = Vs'*Sb*Vs;
%    
%   
%    % Diagonalize the spin blocked hamiltonian
%    [Vbb, Dbb] = eig(Hbb);
%    
%    
%    if number==3
%        figure(1)
%        colormap hot
%        subplot(1,2,1)
%        imagesc(Hbb)
%        subplot(1,2,2)
%        imagesc(abs(Vbb))
%        s1 = zeros(540,1); s1(50,1) = 1/sqrt(2); s1(51,1) = 1/sqrt(2) 
%        s2 = zeros(540,1); s2(100,1) = 1;
%        s3 = zeros(540,1); s3(200,1) = 1;
%        s4 = zeros(540,1); s4(350,1) = 1;
%        s5 = zeros(540,1); s5(450,1) = 1;
%        s6 = zeros(540,1); s6(510,1) = 1;
%        
%        Vbb(:,1)'*Sbb*Vbb(:,1)
%        s1'*Sbb*s1
%        s2'*Sbb*s2
%        s3'*Sbb*s3
%        s4'*Sbb*s4
%        s5'*Sbb*s5
%        s6'*Sbb*s6
%    end
%    
   
   % Save the eigen vectors and eigen values of the Hamiltonian
%    V_tot(posi:posf,posi:posf) = Vbb;
%    D_tot(posi:posf,posi:posf) = Dbb;
   
end


end

