clc; clear all; close all;
% Get pauli matrices
sigx = pauli(1) ; sigy = pauli(2) ; sigz = pauli(3) ; 

% Obtain 1D Heisenberg Hamiltonian with periodic boundary condition
Nx = 12 ; Ny = 1 ; 
Jx = 1 ; Jy = Jx ; Jz = Jx ; 
pbc = 1 ;  % Periodic boundary condition
Hh = H_heisenbergp(Nx, Ny, Nx*Ny, Jx, Jy, Jz, pbc) ; 
Ntot = Nx*Ny ; 

% Choose J positive so ground state is antiferromagnet
Htot = Hh/4 ;  

[V, E] = eig(full(Htot)) ; 

Egs = diag(E) ; 
Egs(1:6)
% Check if eigenstate is normalized
% V(:,1)'*V(:,1)  

% Take |g.s.> and construct \rho
rho_gs = V(:,1)*V(:,1)' ; 

% Get rho for last spins
% Obtain rho of first and last spins
rho_s = rho_gs ; 
for middle_s = 1:Nx*Ny-2
  rho_m = 0 ;
  rho_m = trace_middle_spin(rho_s) ;
  rho_s = rho_m ; 
end
rho_1_l =  rho_s ;

% Get rho for 2nd and 3rd spin
Ntot = Nx*Ny ; 

% Array to save Sx, Sy, Sz for each spin
% rows are spin number and columns are Sx, Sy, Sz

% Average spin components
Spin = zeros(Ntot, 3) ; 

% Entropy of each local spin
S_en = zeros(Ntot, 1) ; 

for sp  = 1: Ntot
  rho_sp = 0 ; 
  rho_sp = rho_loc(rho_gs, sp, Ntot) ; 
  Spin(sp, 1) = trace(rho_sp*sigx) ; 
  Spin(sp, 2) = trace(rho_sp*sigy) ; 
  Spin(sp, 3) = trace(rho_sp*sigz) ; 
  S_en(sp, 1) = get_entropy(rho_sp) ; 
end

Spin 

% Compute mutual information between 1 and 4 th spins
S_1_l = get_entropy(rho_1_l) ; 
M_1_l = S_en(1,1) + S_en(Ntot, 1) - S_1_l ;

% Mutual information is nonzero between first and last spins 
% To verify if code is computed correctly, I traced last spin to get rho1 and
% also traced first spin to get rhoN and they give correct result. 
dlmwrite('12Spin.txt', Spin) ; 
dlmwrite('12Spin_entropy.txt', S_en) ;


% Compute entropy for half chain
% Trace out half spin chain dof
rho_hlf = rho_half(rho_gs, Ntot) ; 
 
% Compute entropy 
S_hlf = get_entropy(rho_hlf) 


return
% The latter code needs to run on kohn

% Compute overlap with Valencebond state
% Singlet state :
up = [1 ; 0] ; dwn = [0; 1] ; 
snglt_st = (kron(up, dwn) - kron(dwn, up))/sqrt(2) ; 

for ii = 1:(Ntot/2-1)
    snglt_st = sparse(kron(snglt_st, snglt_st)) ;
end 
% 
amp_VBS = snglt_st'*V(:,1) ; 
proj_VBS = amp_VBS'*amp_VBS

E_VBS = snglt_st'*Htot*snglt_st
