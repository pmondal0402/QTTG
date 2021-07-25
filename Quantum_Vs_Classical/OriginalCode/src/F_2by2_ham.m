function H = F_2by2_ham(Jh,Nspins)
dim = 2^Nspins;
%pauli matrices
sigmax = [0 1;1 0]/2;
sigmay = [0 -1i;1i 0]/2;
sigmaz = [1 0;0 -1]/2;

% pair pauli matrices
Sx = kron(sigmax,sigmax);
Sy = kron(sigmay,sigmay);
Sz = kron(sigmaz,sigmaz);

% hamiltonian construction
Hx = zeros(dim,dim);
Hy = zeros(dim,dim);
Hz = zeros(dim,dim);
 
for a=1:Nspins-1
   nl = a-1;
   nr = Nspins - a -1;
   diml = 2^nl;
   dimr = 2^nr;
   
   Hx = Hx + kron(kron(eye(diml),Sx),eye(dimr));
   Hy = Hy + kron(kron(eye(diml),Sy),eye(dimr)); 
   Hz = Hz + kron(kron(eye(diml),Sz),eye(dimr)); 

end
H = -Jh*(Hx + Hy + Hz);
end
