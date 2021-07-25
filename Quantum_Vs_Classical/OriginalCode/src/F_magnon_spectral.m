function [states_plus, val] = F_magnon_spectral(Egs,GST,V,D,Nspins,N,n,s)

% Hilbert Space Dimensions 
dim_e = 4^N;
dim_sp = n^Nspins;
dim = dim_e*dim_sp;

val = 0;

% Call a-operator
a = F_Aop(n);

% Construct the ladder operators
af = zeros(dim_sp,dim_sp,Nspins);
adagf = zeros(dim_sp,dim_sp,Nspins);

for c=1:Nspins
    nl = (c-1);
    nr = Nspins - c;
    diml = n^nl;
    dimr = n^nr;

    af(:,:,c) = kron(kron(eye(diml),a),eye(dimr));
    adagf(:,:,c) = kron(kron(eye(diml),a'),eye(dimr));
end

% Evaluate the Weight of each state 
c1 = 1;
%c2 = 1;
for k=1:length(GST)
   % Pre-factors of delta function for each states 
    wkp = 0;
    wkm = 0;
   
    for a=1:N
     % Magnonic creation and annhilation operators
     adag_op = kron(eye(dim_e),adagf(:,:,a));
     a_op = kron(eye(dim_e),af(:,:,a)); 
     
     % Operators acting on the ground state
     adagGS = adag_op*GST;
     aGS = a_op*GST;
     
     %Amplitude with Psi-K vectors 
     wkc = (V(:,k))'*adagGS;
     wka = (V(:,k))'*aGS;
     
     % Evaluate Weights for a k-mode
     weight_1 = wkc'*wkc;   
     weight_2 = wka'*wka;
     
     % Add for the total mode
     wkp = wkp + weight_1;   
     wkm = wkm + weight_2;  

    end
   
   val = val + wkp - wkm;
    % If the weight is non-zero then note down the state number, its energy and the weight
   if wkp>10^(-3)
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
   
   %if wkm>10^(-3)
   %  state_no = k;
   %  energy = D(k,k);
   %  peak_pos = -1*(energy-Egs);
   %  weight = wkm;

   %    states_minus(c2,1) = state_no;
   %    states_minus(c2,2) = energy;
   %    states_minus(c2,3) = peak_pos;
   %    states_minus(c2,4) = weight;     
   %   c2 = c2 + 1;
   %end

end

end
