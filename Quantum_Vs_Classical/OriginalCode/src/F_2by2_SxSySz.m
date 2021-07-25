function Sa = F_2by2_SxSySz(site,Nspins,choice)
dim = 2^Nspins;

%pauli matrices
sigmax = [0 1;1 0]/2;
sigmay = [0 -1i;1i 0]/2;
sigmaz = [1 0;0 -1]/2;

% All matrices put together
sall = zeros(2,2,3);
sall(:,:,1) = sigmax;
sall(:,:,2) = sigmay;
sall(:,:,3) = sigmaz;

% hamiltonian construction
Sa = zeros(dim,dim);

   nl = site-1;
   nr = Nspins - site;
   diml = 2^nl;
   dimr = 2^nr;

   Sa = kron(kron(eye(diml),sall(:,:,choice)),eye(dimr));

end

