function [states_plus, states_minus, val] = F_spectral(Egs,GST,V,D,Nspins,N,n,s)

% Hilbert Space Dimensions 
dim_e = 4^N;
dim_sp = n^Nspins;
dim = dim_e*dim_sp;


val = 0;

% Evaluate the Weight of each state 
c1 = 1;
c2 = 1;
for k=1:length(GST)

    wkp= 0;
    wkm= 0;
    for a=1:N

           % Electronic-creation annhilation operators
           Cdag_up = sparse(kron(F_Cdag_up(a,N),eye(dim_sp)));
           Cdag_dn = sparse(kron(F_Cdag_dn(a,N),eye(dim_sp)));
           C_up = Cdag_up';
           C_dn = Cdag_dn';
           
	       % Operators acting on the ground state
           Gpu = Cdag_up*GST;
           Gpd = Cdag_dn*GST;
           Gmu = C_up*GST;
           Gmd = C_dn*GST;

           % Amplitude with Psi_k states (of up and down electrons)
           Cpu = (V(:,k))'*Gpu; Cpd=(V(:,k))'*Gpd;
           Cmu = (V(:,k))'*Gmu; Cmd = (V(:,k))'*Gmd;

           % Now evaluate the weights for a k-mode for this site
           wpu = Cpu'*Cpu;
           wmu = Cmu'*Cmu;
           wpd = Cpd'*Cpd;
           wmd = Cmd'*Cmd;

           % Add the weights for a k-mode for all the sites
           wkp = wkp + wpu + wpd; %( For plus peaks)
           wkm = wkm + wmu + wmd; %( For minus peaks) 	   
	   
    end
   
   val = val + wkp + wkm;
   
   % If the weight is non-zero then note down the state number, its energy and the weight
   if wkp>10^(-4)
     state_no = k;
     energy = D(k,k);
     peak_pos = energy-Egs;
     weight = wkp;

       states_plus(c1,1) = state_no;
       states_plus(c1,2) = energy;
       states_plus(c1,3) = peak_pos;
       states_plus(c1,4) = weight; 
       c1 = c1+1;
   end
   
   if wkm>10^(-4)
     state_no = k;
     energy = D(k,k);
     peak_pos = -1*(energy-Egs);
     weight = wkm;

       states_minus(c2,1) = state_no;
       states_minus(c2,2) = energy;
       states_minus(c2,3) = peak_pos;
       states_minus(c2,4) = weight;     
       c2 = c2 + 1;
   end

end

end
