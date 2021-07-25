function avg =  F_aadag_avg_spin(psi_t,Nspins,n,s,nterms)
% Call a-operator
a = F_Aop(n);

% Construct the ladder operators
dim = n^Nspins;
af = zeros(dim,dim,Nspins);
adagf = zeros(dim,dim,Nspins);

for c=1:Nspins
    nl = (c-1);
    nr = Nspins - c;
    diml = n^nl;
    dimr = n^nr;

    af(:,:,c) = kron(kron(eye(diml),a),eye(dimr));
    adagf(:,:,c) = kron(kron(eye(diml),a'),eye(dimr));
end

% Construct the Hamiltonian
avg = zeros(1,3*Nspins);
for c=1:Nspins
   
   Spi = sqrt(2*s)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s)*af(:,:,c);
   Smi = sqrt(2*s)*adagf(:,:,c)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s);

   Sxi =  (Spi + Smi)/2;
   Syi =  (Spi - Smi)/(2*1i);
   Szi = (s*eye(dim) - adagf(:,:,c)*af(:,:,c));

   pos = (c-1)*3 + 1;
   avg(1,pos) = psi_t'*Sxi*psi_t;
   avg(1,pos+1) = psi_t'*Syi*psi_t;
   avg(1,pos+2) = psi_t'*Szi*psi_t;

end

end
