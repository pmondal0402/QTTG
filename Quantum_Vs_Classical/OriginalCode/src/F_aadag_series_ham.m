function H = F_aadag_series_ham(Jh,Nspins,n,s,nterms,ni)

ani = 0*10^(-3);

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
H = sparse(dim,dim);
for c=1:Nspins-1

Spi = sqrt(2*s)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s)*af(:,:,c);
Smi = sqrt(2*s)*adagf(:,:,c)*F_sqrt_series(af(:,:,c),adagf(:,:,c),nterms,s);
Szi = (s*eye(dim) - adagf(:,:,c)*af(:,:,c));

Spii = sqrt(2*s)*F_sqrt_series(af(:,:,c+1),adagf(:,:,c+1),nterms,s)*af(:,:,c+1);
Smii = sqrt(2*s)*adagf(:,:,c+1)*F_sqrt_series(af(:,:,c+1),adagf(:,:,c+1),nterms,s);
Szii = (s*eye(dim) - adagf(:,:,c+1)*af(:,:,c+1));

if ni==0
    H = H - Jh*0.5*(Spi*Smii + Smi*Spii) - (Jh+ani)*Szi*Szii;
else
    H = H - Jh*0.5*(Spi*Smii + Smi*Spii) - (Jh+ani)*Szi*Szii + (Jh+ani)*adagf(:,:,c)*af(:,:,c)*adagf(:,:,c+1)*af(:,:,c+1);
end
end
H = sparse(H);
end
