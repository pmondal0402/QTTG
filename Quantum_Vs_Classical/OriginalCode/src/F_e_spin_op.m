function Sop = F_e_spin_op(site,C_up,C_dn,N,choice)
% electronic single particle noninteracting Hamiltoinan
H = zeros(N,N);
for a=1:N-1
   H(a,a+1) = -1;
   H(a+1,a) = -1;
end

[V, D] = eig(H); 

% Construct the U-operators
U = zeros(2,2,N,N);

for n=1:N
   for np=1:N
        U(1,1,n,np) = V(n,np);
        U(2,2,n,np) = V(n,np);
   end
end
 

% Construct the pauli matrices
sigmax = 0.5*[0 1;1 0];
sigmay = 0.5*[0 -1i;1i 0];
sigmaz = 0.5*[1 0;0 -1];

sigma = zeros(2,2,3);
sigma(:,:,1) = sigmax;
sigma(:,:,2) = sigmay;
sigma(:,:,3) = sigmaz;

% Now construct the spin operators
dim = length(C_up);
Sop = zeros(dim,dim);

for n=1:N
    for np=1:N
       mat = (U(:,:,site,np))'*sigma(:,:,choice)*(U(:,:,site,n));
       Sop = Sop + mat(1,1)*(C_up(:,:,np))'*(C_up(:,:,n)) ...
                 + mat(1,2)*(C_up(:,:,np))'*(C_dn(:,:,n)) ...
                 + mat(2,1)*(C_dn(:,:,np))'*(C_up(:,:,n)) ...
                 + mat(2,2)*(C_dn(:,:,np))'*(C_dn(:,:,n));
    end
end

end
