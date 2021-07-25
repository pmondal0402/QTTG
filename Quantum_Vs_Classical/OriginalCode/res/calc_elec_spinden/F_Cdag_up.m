function op = F_Cdag_up(a,N)
% Dimension of Hilbert Space
dim = 4^N;
Cdag = [0 0 0 0;...
      1 0 0 0;...
      0 0 0 0;...
      0 0 1 0];

P = [1 0 0 0; 
     0 -1 0 0; 
     0 0 -1 0; 
     0 0  0 1];

% Construct the permutation operator at point `a'
run = a-1;
perm = 1;
for c=1:run
    perm = kron(perm,P);
end

% Construct the operator at point `a'
op = kron(perm,Cdag);
dim_r = 4^(N-a);
op = kron(op,eye(dim_r));

end
