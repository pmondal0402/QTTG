clc
clear all

tag = 4;
%------------
% Variables
%------------
% Hisenberg Coupling
Jh = 1;
% Sd-coupling
Jsd = 1;
% Total number of spins
Nspins = 2;
% Trunction of HP-Boson
n = 3;
% Spin Value of HP-Boson
s = 5/2;
% Total number of sites
N = Nspins;
% Hopping potential
gamma = 1;
% Initial time
ti = 0;
% final time
tf = 5;
% time step
dt = 0.1;
% Plancks Constant
hbar = 0.658211951;

%---------------------------------------------
% Inital confuguration of spins and electrons
%---------------------------------------------
% Angles of all local spins
theta_sp = [pi 0 pi/2];
phi_sp = [0 0 0];

% State of electrons (1=vacc,2=up,3=dn,4=updn) 
state_e = [2 1 3];

% Dimensions
dim_e = 4^N;
dim_sp = n^Nspins;

%---------------
% Initial State
%---------------
% State vector initiation of local spins
psi_sp = F_aadag_initial(theta_sp,phi_sp,Nspins,n,s);
% State vector initiation of electrons
psi_store = zeros(4,N);
psi_e = 1;

for a=1:N
    pos = state_e(1,a);
    psi_store(pos,a) = 1;
end

for a=1:N
    psi_e = kron(psi_e,psi_store(:,a)); 
end

% Total State vector initialization
psi = kron(psi_e,psi_sp);

%--------------
% Hamiltonian
%--------------
% Spin only Hamiltonian
H_sp = F_aadag_ham(Jh,Nspins,n,s); H_sp = kron(eye(dim_e),H_sp);
% Electronic Hamiltonian
H_e =  F_ccdag_ham(gamma,N); H_e = kron(H_e,eye(dim_sp));
% Sd-interation hamiltonian
H_jsd = F_cacadag_ham(Jsd,Nspins,N,n,s);
% Total Hamiltonian
H = H_sp + H_e + H_jsd;

%----------------
% Time Evolution
%----------------
% Average local spins
avg_sp = zeros(length(ti:dt:tf),3*Nspins);
% Average electron spin density
avg_e_spd = zeros(length(ti:dt:tf),3*Nspins);
% Average electron charge density
avg_e_ch = zeros(length(ti:dt:tf),Nspins);
% Average Magnon number
avg_mag_no = zeros(length(ti:dt:tf),1);

c=1;
for t=ti:dt:tf
   t
   psi_t = expm(-1i*H*t/hbar)*psi;
   [avg0 magno] = F_fock_avg_sp(psi_t,N,Nspins,n,s);
   [avg1 avg2] = F_fock_avg_e(psi_t,N,Nspins,n,s);
   avg_sp(c,:) = real(avg0);
   avg_mag_no(c,1) = magno;
   avg_e_spd(c,:) = real(avg1);
   avg_e_ch(c,:) = real(avg2);
   c=c+1;
end

%-------------
% Save data 
%-------------
vars = [N Jh Jsd];

string = ['../../results/test/emagnon_hp/vars' num2str(tag) '.dat'];
F_print_m(string,vars,4);

string = ['../../results/test/emagnon_hp/avg_sp' num2str(tag) '.dat'];
F_print_m(string,[(ti:dt:tf)', avg_sp],4);

string = ['../../results/test/emagnon_hp/avg_e_spd' num2str(tag) '.dat'];
F_print_m(string,[(ti:dt:tf)', avg_e_spd],4);

string = ['../../results/test/emagnon_hp/avg_e_ch' num2str(tag) '.dat'];
F_print_m(string,[(ti:dt:tf)', avg_e_ch],4);

string = ['../../results/test/emagnon_hp/avg_mag' num2str(tag) '.dat'];
F_print_m(string,[(ti:dt:tf)', avg_mag_no],4);

