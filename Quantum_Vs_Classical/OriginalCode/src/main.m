% This test is for 1D AFM 
% Checking if GS is ordered for high spin value ?
% Ans : 

clc
clear all
close all ; 
% Controlling number of cores to use
maxNumCompThreads(25);
tag = 4;
%------------
% Variables
%------------
% Hopping potential
gamma = 1;
% Hisenberg Coupling
% -ve for AFM; +ve for FM
Jh = -0.01;
% Sd-coupling
Jsd = 0.01 ;

% Total number of spins
Lx = 4 ; 
Ly = 1 ; 
Nspins = Lx*Ly ;

% State of electrons (1=vacc,2=up,3=dn,4=updn) 
state_e = [2 1 1 1];

% Total number of sites
N = Nspins ;
sv_frm = 0 ; 

% TODO : Change spin value below 
for val = 1:1 % 0.5:8% 0.5:0.5:3 
   sv_frm = sv_frm  + 1 ; 


   % Spin Value of HP-Boson
   s = val 
   % Trunction of HP-Boson
   n = 2*s+1 ;
   Bz = -0.06  ; % TODO Set val to -0.05
   B0x = -0.01 ; % TODO Set val to -0.01
   Tw = 15 ; t0 = 50 ; 
   % Plancks Constant
   hbar = 0.658211951;

   %---------------------------------------------
   % Inital confuguration of spins and electrons
   %---------------------------------------------

   % Angles of all local spins
   theta_sp = [pi 0 pi/2];
   phi_sp = [0 0 0];
   
   
   % Dimensions
   dim_e = 4^N;
   dim_sp = n^Nspins;
   
   %---------------
   % Initial State
   %---------------
   [Jx,Jy,Jz,Jminus,Jplus] = Joperators(s) ; 
   for nval = 1:n
     Jz(nval, nval) = s - (nval-1) ;
   end

   [V, E] = eig(Jx) ; 
   [Vz, Ez] = eig(Jz) ; 
  
   J1x = kron(Jx, speye(n^(Nspins-1) ) ) ;  
   J1y = kron(Jy, speye(n^(Nspins-1) ) ) ;  
   J1z = kron(Jz, speye(n^(Nspins-1) ) ) ;  

   J2x = kron( kron(speye(n), Jx), speye(n^(Nspins-2) ) ) ;  
   J2y = kron( kron(speye(n), Jy), speye(n^(Nspins-2) ) ) ;  
   J2z = kron( kron(speye(n), Jz), speye(n^(Nspins-2) ) ) ;  

   J3x = kron( kron(speye(n^2), Jx), speye(n^(Nspins-3) ) ) ;  
   J3y = kron( kron(speye(n^2), Jy), speye(n^(Nspins-3) ) ) ;  
   J3z = kron( kron(speye(n^2), Jz), speye(n^(Nspins-3) ) ) ;  

   J4x = kron( kron(speye(n^3), Jx), speye(n^(Nspins-4) ) ) ;  
   J4y = kron( kron(speye(n^3), Jy), speye(n^(Nspins-4) ) ) ;  
   J4z = kron( kron(speye(n^3), Jz), speye(n^(Nspins-4) ) ) ;  

   psi_sp = V(:, end) ; 
   psi_spz = Vz(:, 1) ; % State of spin down
   psi_spz2 = Vz(:, end) ; % State of spin up
    
   % Total State vector initialization
   % TODO : Compute the state manually to avoid negative sign issue

 
   Neel_vec = kron(psi_spz2, kron(psi_spz, kron(psi_spz2, psi_spz))) ; 
   Neel_v2 = (kron(psi_spz, psi_spz2) - kron(psi_spz2, psi_spz))/sqrt(2) ; 
   proj_neel = Neel_vec*sparse(Neel_vec') ; 
  
   %--------------
   % Hamiltonian
   %--------------
   % Spin only Hamiltonian
   [H_sp, S1x] = F_aadag_ham(Jh,Nspins,n,s, Bz, Lx, Ly);  

   % Utkarsh's original code is used to make sure I didn't add any bug
   % [H_sp] = F_aadag_ham(Jh,Nspins,n,s);  
   % Ans. It is same.

   % Add B0x field to S1
   Hpulse = B0x*S1x ;

   % Total Hamiltonian
   H = H_sp + Hpulse;
   
   % Is the |GS> unique in FOCK space for even number of spins ?
   % Yes ! It is unique for all spin values 
   if val < 2
     k = length(H) ; 
   else
     k = 100 ; 
   end

   [V_exact, E_exact] = eigs(H, k) ; 

   % [V_exact, E_exact] = eig(H) ; 
   E_exact_val = diag(E_exact) ;
 
   % Compute local spin polarization 
   rho_gs = V_exact(:, 1)*sparse(V_exact(:, 1)') ;  

   S1(sv_frm, 1) = trace(rho_gs*J1x)/s ; 
   S1(sv_frm, 2) = trace(rho_gs*J1y)/s ; 
   S1(sv_frm, 3) = trace(rho_gs*J1z)/s ; 

   S2(sv_frm, 1) = trace(rho_gs*J2x)/s ; 
   S2(sv_frm, 2) = trace(rho_gs*J2y)/s ; 
   S2(sv_frm, 3) = trace(rho_gs*J2z)/s ; 


   S3(sv_frm, 1) = trace(rho_gs*J3x)/s ; 
   S3(sv_frm, 2) = trace(rho_gs*J3y)/s ; 
   S3(sv_frm, 3) = trace(rho_gs*J3z)/s ; 

   S4(sv_frm, 1) = trace(rho_gs*J4x)/s ; 
   S4(sv_frm, 2) = trace(rho_gs*J4y)/s ; 
   S4(sv_frm, 3) = trace(rho_gs*J4z)/s ; 
   % Save rho 
   % save(sprintf('../res/saved_rho/rho_gs%03d.mat', sv_frm),'rho_gs');
   % save(sprintf('../res/saved_rho/eigval%03d.mat', sv_frm),'E_exact_val');
end

S1
S2
S3
S4

trace(rho_gs*proj_neel)
% Save rho to check contribution from Neel state
% save('../res/saved_rho_gs.mat', 'rho_gs') ; 

% Electronic Hamiltonian
H_e =  F_ccdag_ham(gamma,N); H_e = kron(H_e,speye(dim_sp));

% Sd-interation hamiltonian
H_jsd = F_cacadag_ham(Jsd,Nspins,N,n,s);

% Total Hamiltonian modified (see line 112) 
H_tot = kron(speye(dim_e), H) + H_e + H_jsd ; 

% TODO : Add psi_e to psi evol

% State vector initiation of electrons
psi_store = sparse(4,N);
psi_e = 1;

for a=1:N
    pos = state_e(1,a);
    psi_store(pos,a) = 1;
end

for a=1:N
    psi_e = kron(psi_e,psi_store(:,a)); 
end

% Time evolution
% Initial time
ti = 0 ; 
% final time
tf = 200 ; 
% time step
dt = 0.01 ; 
% File to save data
file_sp = fopen('../res/spinQM.txt','w') ;
% Setting initial state psi
psi = kron(psi_e, Neel_vec) ;

% Modify J1x, J1y to include elec dof  ...
J1x = kron(speye(dim_e), J1x) ; J1y = kron(speye(dim_e), J1y);
J1z = kron(speye(dim_e), J1z) ; 

J2x = kron(speye(dim_e), J2x) ; J2y = kron(speye(dim_e), J2y);
J2z = kron(speye(dim_e), J2z) ; 

J3x = kron(speye(dim_e), J3x) ; J3y = kron(speye(dim_e), J3y);
J3z = kron(speye(dim_e), J3z) ; 

J4x = kron(speye(dim_e), J4x) ; J4y = kron(speye(dim_e), J4y);
J4z = kron(speye(dim_e), J4z) ;
 
c = 1 ;  
for t=ti:dt:tf
  t
  [psi] = get_CN(sparse(H_tot), sparse(psi), dt, 0, 1) ;

   %{ 
   [avg0 magno] = F_fock_avg_sp(psi,N,Nspins,n,s);
   [avg1 avg2] = F_fock_avg_e(psi,N,Nspins,n,s);
   avg_sp(c,:) = real(avg0);
   avg_mag_no(c,1) = magno;
   avg_e_spd(c,:) = real(avg1);
   avg_e_ch(c,:) = real(avg2);
   %}

 % TODO : replace last 6 operators in argument of function below to get 
 % corresponding observable 
  avg_sp1(c, :) = measr_spin(psi, J1x, J1y, J1z,...
                                 J2x, J2y, J2z,...
                                 J3x, J3y, J3z,...
                                 J4x, J4y, J4z) ; 
  % Save data
  As = [t; real(avg_sp1(c, :))'] ; 
  fprintf(file_sp,'%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n', As) ; 

  % Save psi
  save(sprintf('../res/saved_rho/psi_t%03d.mat', c),'psi');
  % Add + 3 for next spin
  testSx(c) = avg_sp1(c, 1+3) ; 
  testSy(c) = avg_sp1(c, 2+3) ; 
  testSz(c) = avg_sp1(c, 3+3) ; 
  %}
  time(c) = t ; 
  c = c + 1 ; 
end

return

% Is it entangled ? 
plot(time, testSx/s, 'b-')
hold on
plot(time, testSy/s, 'r-')
plot(time, testSz/s, 'k-')
plot(time, sqrt(testSx.^2 + testSy.^2 + testSz.^2)/s)
hold off

