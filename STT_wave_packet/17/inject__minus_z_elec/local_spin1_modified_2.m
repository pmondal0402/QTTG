%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define initial parameters and matrices
%
% Set solver
% direct solver ( i.e. \)
dirsolver = 0;

% Iterative solver
itersolver = 1;

% Spectral decomposition 
spectral_decomp = 0;

Purity_state = 1;

%Decides which measurements to take

Partial_trace_conf_method = 1;

loc_purity = 1;			% Purity of local spins

electron_state = 1;
loc_spin1_state = 1;
loc_spin_middle = 1;
loc_spin_last_state = 0;
movie = 0;




% #############################################################################
% Not using the package here 
%addpath(genpath('/mnt/graphene/home/pmondal/QM_transport/project5/step_3/Purity/QETLAB-0.9'))
%addpath(genpath('/home/pmondal/test/Purity/Purity/Purity/QETLAB-0.9'))
%addpath(genpath('/home/pmondal/test/local_spin/QETLAB-0.9'))
%addpath((genpath('/mnt/graphene/home/pmondal/QM_transport/project5/step_3/Purity/QETLAB-0.9')))

%% #############################################################################


% Declaring initial parameters i.e. length of chain, wavepacket parameters  etc.............

Lx = 401;
Ly = 1;
N_loc_spin = 5;     % No. of local spins in System... 
pos_loc_i = 198;%25*Ly+1;%1;%10*Ly+1;%15*Ly+1;%24*Ly + 1;%65*Ly + 1; %13*Ly + 1;%70*Ly + 1 ;      % Position where local spin starts

dt = 0.01; % Iterative solver time step
tmax = 90;
% Tight binding and rashba hopping rates .........

t = 1;

t_so = 0.;%0.1;
lx1 = 0;%11;        % where Rashbha starts
lx2 = 4;%41;        % where Rashbha ends

a = 1 ; %not sure about this
%kx = 1.5;%pi/(2*a);%2.6/a;
kx =  0.1;%0.5;% 0.105;%0.0055;
%deltakx = 2./a;
deltakx = 0.2;%0.2;%0.0805;
Jsd =  -0.1;%-0.5;%*t;%2;%1.0;
Jheis = -0.1;%*t;%0.1;%0.1;
Jx = Jheis;
Jy = Jheis;
Jz = 1.005*Jheis;
Bz = -0;
lambda = 0.0;
%Bz = -5.1;

theta = 0;
phi = 0;
LU_decomp = 0;




sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

sigmax = [0 1;1 0];
sigmay = [0 -1i;1i 0];
sigmaz = [1 0;0 -1];
spin_conf = 2^(N_loc_spin + 1);

rashbha_x = (-t*eye(2) -1i*t_so*sigma_y);
rashbha_y = (-t*eye(2) +1i*t_so*sigma_x);


%fileID2 = fopen('site_vs_Py.txt','w');
%fileID3 = fopen('site_vs_Pz.txt','w');
%fileID4 = fopen('site_vs_P.txt','w');

%file = fopen('data_for_cpu_time_Lx_Ly_N_loc_spin_size(Hamilto)_nnz(Hamilto)_cput_cpu2_time_cpu_trace_ovr_spin.txt','w');
 fileID_e = fopen('electron_data.txt','w');

%fileID1 = fopen('P1_all.txt','w');
%fileID2 = fopen('P2_all.txt','w');
%fileID3 = fopen('P3_all.txt','w');
%fileID4 = fopen('P4_all.txt','w');
%fileID5 = fopen('P5_all.txt','w');
%fileID6 = fopen('P6_all.txt','w');
%fileID7 = fopen('P7_all.txt','w');
%fileID8 = fopen('P8_all.txt','w');
%fileID9 = fopen('P9_all.txt','w');
%fileID10 = fopen('P10_all.txt','w');

% 
%  % spin = 1;
   fileID1h = fopen('Purity_2nd_loc.txt','w');
%   fileID1h = fileIDh;
fileID2h = fopen('Purity_3rd_loc.txt','w');
%fileID3h = fopen('P3_h.txt','w');
%fileID4h = fopen('P4_h.txt','w');
%fileID5h = fopen('P5_h.txt','w');
%fileID6h = fopen('P6_h.txt','w');
%fileID7h = fopen('P7_h.txt','w');
%fileID8h = fopen('P8_h.txt','w');
%fileID9h = fopen('P9_h.txt','w');
%fileID10h = fopen('P10_h.txt','w');



fileID1 = fopen('Purity_1st_loc.txt','w');
fileID1_2 = fopen('Purity_2nd_loc.txt','w');
fileID1_3 = fopen('Purity_3rd_loc.txt','w');
fileID1_4 = fopen('Purity_4th_loc.txt','w');


%fileID2 = fopen('N_loc_spin_2_tx_Px_Py_Pz_P_sd_no_rashba_plus_Heisenberg_coupling_local_l4_tx=2fs_delta_kx=2.txt','w');
fileID3 = fopen('Purity_last_loc.txt','w');

fileID_p_e = fopen('tx_Purity_elec_sub_space.txt','w');
fileID_p_loc = fopen('tx_Purity_loc_sub_space.txt','w');


%% Tight Binding Hamiltonian 

%% Define position space hamiltonians.... Upper half..TB Hamiltonian
% H1................................................

H1_pos = sparse(Lx*Ly,Lx*Ly);
H1_pos_i = sparse(Lx*Ly,Lx*Ly);
H1_pos_o = sparse(Lx*Ly,Lx*Ly);


for n = 1:1:((Lx-1)*Ly)
       H1_pos(n,n+Ly) = t;
      

    if n > (lx1-1)*Ly && n<(lx2*Ly+1) %((lx2+1)*Ly+1)
        H1_pos_i(n,n+Ly) = H1_pos(n,n+Ly);
        H1_pos_o(n,n+Ly) = 0;
    else
        H1_pos_i(n,n+Ly) = 0;
        H1_pos_o(n,n+Ly) = H1_pos(n,n+Ly);
    end    
     
       
end    


% H0................................................


H0_pos = sparse(Lx*Ly,Lx*Ly);
H0_pos_i = sparse(Lx*Ly,Lx*Ly);
H0_pos_o = sparse(Lx*Ly,Lx*Ly);
%Lx

for n = 1:1:(Lx*Ly-1)

    %Defining H0 inside, out................................
    if mod(n,Ly) == 0
        continue
    else
       
        H0_pos(n,n+1) = t;
    end 
    
    if n > (lx1*Ly) && n < ((lx2+1)*Ly+1)
        H0_pos_i(n,n+1) = H0_pos(n,n+1);
        H0_pos_o(n,n+1) = 0;
    else
        H0_pos_i(n,n+1) = 0;
        H0_pos_o(n,n+1) = H0_pos(n,n+1);
    end   
end   


% No rashba hamiltonian 
% Write TB hamiltonian 


% Defining Rashbha region and shifting Hamiltonian to spin space




H1_tot_i = kron(H1_pos_i,sparse(rashbha_x));
H1_tot_o = kron(H1_pos_o,sparse((-t*eye(2))));
H1 = H1_tot_i + H1_tot_o;


H0_tot_i = kron(H0_pos_i,sparse(rashbha_y));
H0_tot_o = kron(H0_pos_o,sparse(-t*eye(2)));
H0 = H0_tot_i + H0_tot_o;

H_tot = H0+ H1;

HH = H_tot + ctranspose(H_tot);
H_tot_rashba_spin = kron(HH,eye(2^N_loc_spin)); 



%% Having classical spins 

poss = sparse(Lx*Ly,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% Building |Psi>

Psi_ns = zeros(Lx*Ly,1);

%txx = zeros(tmax,1);
index=1;
for x=1:Lx
    for y=1:Ly
        Psi_ns(index,1) = exp(1i*kx*x - deltakx^2*(x-160)^2/4);
        index=index+1;
    end
end

Psi_ns = Psi_ns/sqrt(Psi_ns'*Psi_ns);
Psi_e=kron(Psi_ns(:,1),([0 1]')) ;%+  kron(Psi_ns(:,1),[0 1]');     % x_spin polarized electron wave packet
Psi_e = Psi_e/sqrt(Psi_e'*Psi_e);


% -------------------
% PSI INJECTED LATER
% -------------------

Psi_ns1 = zeros(Lx*Ly,1);

%txx = zeros(tmax,1);
index1 = 1;
for x=1:Lx
    for y=1:Ly
        Psi_ns1(index1,1) = exp(1i*kx*x - deltakx^2*(x-180)^2/4);
        index1=index1+1;
    end
end

Psi_ns1 = Psi_ns1/sqrt(Psi_ns1'*Psi_ns1);
Psi_e1=kron(Psi_ns1(:,1),([0 1]')) ;%+  kron(Psi_ns(:,1),[0 1]');     % x_spin polarized electron wave packet
Psi_e1 = Psi_e1/sqrt(Psi_e1'*Psi_e1);







%sparse(Psi_e*Psi_e')
if N_loc_spin ~= 0
    %Psi = kron( kron(Psi_e,[1 1]'), [ 1 1]');    % two local spins... Up spin
    %local electrons 
    Psi = kron(Psi_e, spinor1(theta,phi));
    for i =1: N_loc_spin-1                     % N_loc_spin no. of local spins
        Psi = kron(Psi, spinor1(theta,phi));
    end
%{  
  Psi=Psi/sqrt(Psi'*Psi);

      %  Psi = kron(Psi_e, [1 0]');
    for i = 5:8% N_loc_spin-1                     % N_loc_spin no. of local spins
        Psi = kron(Psi, [1 0]');
    end
%}
    Psi=Psi/sqrt(Psi'*Psi);

else
    Psi=Psi_e/sqrt(Psi_e'*Psi_e);
end

%{
% ------------------
%ABSORBING POTENTIAL
% ------------------
c1 = 80;	%Absorbing potential centered at x = c
a1 = (1 - 16/c^3);
b1 = (1 - 17/c^3)*(1/c^2);
r11 = 75
r22 = 85
width = 10;		% r22 - r11 = 0
del_k_min = c1/(2*width);
E_minim = 0.01;

for sitee = 1:Lx
	xx(sitee) = 2*del_k_min*(sitee - r11);
	yy(sitee) = a*xx(sitee) - b*(xx(sitee))^3 + 4/(c1 - xx(sitee))^2 ...
			-4/(c1 + xx(sitee))^2;
end
%}
% ------------------

%%########################################################################

% Total Hamiltonian 


%h = H_tot_rashba_spin -1i*kron( absorb_pot(100,10,75,85,55),eye(2)); 
  % gaussian = exp(-0.5*(t-t0)^2/Tw^2);
h=  kron(eye(2*Lx*Ly),(H_heisenberg_gen_spin(N_loc_spin,1,N_loc_spin,Jx,Jy,Jz))) ...
    + H_tot_rashba_spin...
    +Jsd*H_sd(Lx,Ly,N_loc_spin,pos_loc_i);% ...
%    + Bz*kron(speye(2*Lx*Ly), applied_field(N_loc_spin)) ;%+Jsd*H_sd_class ;%... 
       % - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sparse(sigma_x),eye(2^(N_loc_spin-1)))))...
%h = kron(h,speye(length(h))) + kron(speye(length(h)), h);
   E_avg =  Psi'*h*Psi
  


   % --------------------
   % local_spin positions
   % --------------------
   local_spin_pos = (pos_loc_i:pos_loc_i + N_loc_spin -1); 
   yy(1:N_loc_spin) = 0;
    
%% .......................................................................
%h = Jsd*H_sd(Lx,Ly,N_loc_spin,pos_loc_i) + H_tot_rashba_spin  +Jheis*kron(eye(2*Lx*Ly),H_heisenberg_gen_spin(N_loc_spin,1,N_loc_spin));
%H = -0.1*kron(eye(2*Lx*Ly),H_heisenberg_gen_spin(5,2,10));
% Psi'*h*Psi
% 
% [V,E] = eigs((h),40,'SA');
% diag(E(1:10,1:10))


 
a= (1+lambda^2);

%{

Hamilto = sparse((1/a)*(h-1i*lambda*(sparse(h)-(Psi'*h*Psi)*(speye(size(h))))));   

Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;

%}



%#############################################################################

t = 0;
it_count = 0;
obs = 0;
%parpool(40)
wave_shoot = 0;

fprintf(1,'Timestepping loop BGN\n');

while t <= 8*tmax		%110*dt --> this time was selected to run time loop only once since obs was taken at 100 timestep 

Hamilto = sparse((1/a)*(h-1i*lambda*(sparse(h)-(Psi'*h*Psi)*(speye(size(h))))));

Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;




 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
		% CPU timing

		cpu1=cputime;

		% Using as solver the direct solver
		if dirsolver ==1
   			Psi = Hp\(Hm*Psi);
		end

		if itersolver ==1
   			brhs = Hm*Psi;
   			[Psi,ierr,relres,niter] = pcg(Hp,brhs);
   			if ierr ~= 0
    				 fprintf(1,'Iterative solver report: ierr = %d \n',ierr);
     				 fprintf(1,'   Timestep:             %d\n',it_count);
     				 fprintf(1,'   Rel residuum:         %e\n',relres);
     				 fprintf(1,'   Number of iterations: %d\n',niter);
     				 pause
   			end
		end

	t = t + dt;
	it_count = it_count + 1;


	cput = cputime - cpu1;

		if it_count==100
  			fprintf(1,'CPU time per time step: %e [s]\n',cput);
		end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if it_count == 100
        wave_shoot = wave_shoot + 1;
    obs = obs + 1;
%    Psi = Psi/sqrt(Psi'*Psi);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


% ############## Trace over configuration space ##########################
%.................Partial trace and density matrices......................

% -------------------------------------------------------------------------------

%Calculating {2^(N_loc_spin+1) x 2^(N_loc_spin+1)} rho matrix in spin space
%i.e. taking trace of configuration space
%rho_tot = sparse(Psi_t*Psi_t');


if N_loc_spin ~= 0

    %% Making Movie 
    if movie == 1
        
        %fprintf('making movie') 
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Psi_movie = zeros(Lx,1);
    
    for j = 1:Lx
        Psi_movie(j,1) = Psi(1+(j-1)*spin_conf: j*spin_conf,1)'...
                        *Psi(1+(j-1)*spin_conf: j*spin_conf,1);
    end
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    sitee = (1:Lx);
    
    ax = plot(sitee,Psi_movie,'b','LineWidth',1.5);
    hold on
    plot(local_spin_pos,yy,'r.','MarkerSize',10)
    hold off
    ylabel('Probability','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')

 ylim([ 0,.2]);
 xlim([0,Lx]);


pbaspect([1.3 1 1]);
lgndc = legend(['t =' num2str(t/0.6591)],'Local Spin','Interpreter','LaTex');
set(lgndc,'color','none','FontName','Times New Roman','FontSize',10)

legend('boxoff')
xlabel('Site','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')
set(ax,'LineWidth',1.5)    
    
    F(obs) = getframe(gcf);
    s = size(F(obs).cdata);
  
        
    end
    
    
    
    
    
    
% Calculating Purity of the system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if Purity_state == 1 

	% Calculating matrix element for |1up><1up|........................................................................



		fprintf('Testing Purity Calculation method \n')


		Purity = sparse(2*Lx*Ly,2*Lx*Ly);
		shift_row = 0;
		for i = 1:2*Lx*Ly

			indx1 = [(1 + shift_row): (2^N_loc_spin + shift_row)];
			shift_colmn = 0;

			for j  = 1:(2*Lx*Ly)

				indx2 = [(1 + shift_colmn) :(2^N_loc_spin + shift_colmn)]; 
				Purity(i,j) = sum(Psi(indx1,1).*conj(Psi(indx2,1)));

				shift_colmn = shift_colmn + 2^(N_loc_spin);
			end
			shift_row = shift_row + 2^N_loc_spin;

		end

		fprintf('Testing Purity of state at a given time\n')
		P = trace(Purity*Purity);
		A_p = [t/0.6591;P;Psi'*Psi;Psi'*h*Psi];
		fprintf(fileID_p_e,'%22.16f %22.16f %22.16f %22.16f\n',A_p);

		% Compare with QETLAB result.......................................................................

	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ..................................................................................................


	rho_spin = sparse(spin_conf,spin_conf);
	cpu2 = cputime;

	if Partial_trace_conf_method == 1


%>>>>>>>>>>>>>>>>> Calculating matrix elements after Partial trace >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	%	fprintf('Partial trace over conf space begins \n')
	%	fprintf('Method 2 \n')
		parfor row = 1:spin_conf


			ixx = [row:spin_conf:size(Psi)];                 %create all indices to your Psi vector (e.g. for given point in config. space)

			for column = 1:spin_conf

				ixx1 = [column:spin_conf:size(Psi)];	% Creating columns for each row ...Here for N_loc_spin no. of spins, no. of columns will be 2^(1+N_loc_spin);
				rho_spin(row,column) = sum(Psi(ixx,1).*conj(Psi(ixx1,1)));

			end

		end

	end

cpu2_time = cputime - cpu2;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Calculating purity of local spin subsystem
if loc_purity == 1
   
    for roww = 1:2^N_loc_spin
        for coll = 1:2^N_loc_spin
                   locc(roww,coll)  = rho_spin(roww,coll) + rho_spin(roww+2^N_loc_spin,coll+2^N_loc_spin);
        end
    end
    
    Purity_locc = trace(locc*locc)
    A_loc = [t/0.6591;Purity_locc];
    
    		fprintf(fileID_p_loc,'%22.16f %22.16f \n',A_loc);
end





% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%fprintf('%s Time taken for Partial trace over configuration space \n',cpu2_time)

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%  Electron spin calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if electron_state == 1

		fprintf('Calculation for electron density begins\n')

% ............................................................................................................................................................
		% Try this code to check block diagnoal elements of identity matrices 
		% Calculation element of |up_e><up_e| ...........................
		rho_electron = sparse(2,2);

		shiftt_elec_row = 0;
		for row_elec = 1:2

			shiftt_elec_colmn = 0;
			for colmn_elec = 1:2

        			 rho_electron(row_elec,colmn_elec) = sum(diag(rho_spin( (1 + shiftt_elec_row) :( 2^N_loc_spin + shiftt_elec_row) ,(1 + shiftt_elec_colmn):(2^N_loc_spin + shiftt_elec_colmn))));
  				 shiftt_elec_colmn = shiftt_elec_colmn + 2^N_loc_spin;

			end

			shiftt_elec_row = shiftt_elec_row + 2^N_loc_spin;

		end

		%fprintf('Testing electron density');

		% Taking observation of electron spin dynamics ............................................................

		 A_e = [t/(0.6591);trace(rho_electron*sigma_x);trace(rho_electron*sigma_y);trace(rho_electron*sigma_z); ...
             ((trace(rho_electron*sigma_x))^2+(trace(rho_electron*sigma_y))^2+(trace(rho_electron*sigma_z))^2)^0.5;trace(rho_electron*rho_electron)];
 		 fprintf(fileID_e,'%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f \n',A_e);



% .........................................................................................................
		full(rho_electron)

% ............................................................................................................................................................

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	% Testing calculation of 1st local spin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if loc_spin1_state == 1

%	fprintf('Starting to calculate spin density of 1st local spin \n');

	% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		% Calculating 1st element of new matrix ..............
 
		shift_rho1_row = 0;
		for row_rho1 = 1:(2^2)

			shift_rho1_colmn = 0;
			for colmn_rho1 = 1:(2^2)	% 2x2 for electron spin and 1st local spin ...Here we are tracing out spins from 2 to N_loc_spin. 

				rho_spin1(row_rho1,colmn_rho1) = sum(diag(rho_spin( (1+ shift_rho1_row) : ( 2^(N_loc_spin-1)+ shift_rho1_row ), (1+ shift_rho1_colmn) :(2^(N_loc_spin -1) + shift_rho1_colmn))));
				shift_rho1_colmn = shift_rho1_colmn + 2^(N_loc_spin -1);
			end

		shift_rho1_row = shift_rho1_row + 2^(N_loc_spin -1);
		end
		% ...................................................



		% .....Tracing electron spin .........................................................
		for l = 1:2
			for m = 1:2
				rho1(l,m) = rho_spin1(l,m) + rho_spin1(l + 2,m+2);
			end
		end

		full(rho1);
		% ...............................................................

		% Taking observation of 1st spin dynamics ...........................................................
		
		A1 = [N_loc_spin/N_loc_spin; t/(0.6591);trace(rho1*sigma_x);trace(rho1*sigma_y);trace(rho1*sigma_z); ...
            ((trace(rho1*sigma_x))^2+(trace(rho1*sigma_y))^2+(trace(rho1*sigma_z))^2)^0.5;trace(rho1*rho1)];
                fprintf(fileID1,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A1);


		% ...............................................................

	% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Calculating inbetween spins .................................................................................................................................

	% Let's start by calculating spin N_loc_spin/2 	

        % Let's start by calculating spin N_loc_spin/2 

 
         if loc_spin_middle == 1

  

		fprintf('Calculating middle spin \n');
for middle_spinn = 2:(N_loc_spin )

	
% Calculating middle spin ................................................

no_sys_aftr_frwd_trce = middle_spinn;%(N_loc_spin/2);



rho_mdl  = zeros(2^no_sys_aftr_frwd_trce,2^no_sys_aftr_frwd_trce);

		for ll1  = 1:length(rho_mdl)
			for mm1 = 1:length(rho_mdl)

				rho_mdl(ll1,mm1) = sum(diag( rho_spin(ll1:length(rho_mdl):length(rho_spin),mm1:length(rho_mdl):length(rho_spin)) ));

			end
		end


% .........................................................................

      rho_middle1 =  zeros(2,2);
      bk_trace_block_size = 2^(no_sys_aftr_frwd_trce - 1);

      rho_middle1(1,1) = sum(diag(rho_mdl(1:bk_trace_block_size,1:bk_trace_block_size ) ));
      rho_middle1(1,2) = sum(diag(rho_mdl(1:bk_trace_block_size, (1+ bk_trace_block_size):2*bk_trace_block_size) ));
      rho_middle1(2,1) = sum(diag(rho_mdl(1 +bk_trace_block_size :2*bk_trace_block_size, 1:bk_trace_block_size) ));

      rho_middle1(2,2) = sum(diag(rho_mdl(1 + bk_trace_block_size:2*bk_trace_block_size, (1+ bk_trace_block_size):2*bk_trace_block_size) ));

    rho_middle1;
             %   A1 = [(N_loc_spin -no_sys_aftr_frwd_trce) ; spin;trace(rho1*sigma_x);trace(rho1*sigma_y);trace(rho1*sigma_z); ((trace(rho1*sigma_x))^2+(trace(rho1*sigma_y))^2+(trace(rho1*sigma_z))^2)^0.5];
              %   fprintf(fileID2,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A1);
 
	% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if middle_spinn == 2
                   P2x = trace(rho_middle1*sigma_x);
                   P2y = trace(rho_middle1*sigma_y);
                   P2z = trace(rho_middle1*sigma_z);
                   P2 = sqrt(P2x^2 + P2y^2 + P2z^2);
                
                  A2 = [middle_spinn;P2x;P2y;P2z;P2];
                  fprintf(fileID1_2,'%10.18f  \t %10.18f \t %10.18f \t %10.18f \t %10.18f\n',A2);



end	

if middle_spinn == 3
                    P3x = trace(rho_middle1*sigma_x);
                    P3y = trace(rho_middle1*sigma_y);
                    P3z = trace(rho_middle1*sigma_z);
                    P3 = sqrt(P3x^2 + P3y^2 + P3z^2);


                  A3 = [middle_spinn;P3x;P3y;P3z;P3];
                  fprintf(fileID1_3,'%10.18f  \t %10.18f \t %10.18f \t %10.18f \t %10.18f \n',A3);

 
end
 
if middle_spinn == 4
                    P4x = trace(rho_middle1*sigma_x);
                    P4y = trace(rho_middle1*sigma_y);
                    P4z = trace(rho_middle1*sigma_z);
                    P4 = sqrt(P4x^2 + P4y^2 + P4z^2);


                  A4 = [middle_spinn;P4x;P4y;P4z;P4];
                  fprintf(fileID1_4,'%10.18f  \t %10.18f \t %10.18f \t %10.18f \t %10.18f\n',A4);


end

if middle_spinn == 5
    rho_last = rho_middle1;
    
    
	% ..................................................................................
	% Taking observation of last spin dynamics ...........................................................		
                      A11 = [N_loc_spin; t/(0.6591);trace(rho_last*sigma_x);trace(rho_last*sigma_y);trace(rho_last*sigma_z); ((trace(rho_last*sigma_x))^2+(trace(rho_last*sigma_y))^2+(trace(rho_last*sigma_z))^2)^0.5;trace(rho_last*rho_last)];
                      fprintf(fileID3,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A11);
  
	% ..................................................................................
%                 
    if wave_shoot == 400
        fprintf('shooting another electron')
        Psi = kron(Psi_e1,spinor(A1(3),A1(4),A1(5),A1(6)));
        Psi = kron(Psi,spinor(A2(2),A2(3),A2(4),A2(5)));
        Psi = kron(Psi,spinor(A3(2),A3(3),A3(4),A3(5)));
        Psi = kron(Psi,spinor(A4(2),A4(3),A4(4),A4(5)));
        Psi = kron(Psi,spinor(A11(3),A11(4),A11(5),A11(6)));
      
   

       % Psi = Psi/sqrt(Psi'*Psi);
        wave_shoot = 0;
    end
        
    
    
    
%                     P4x = trace(rho_middle1*sigma_x);
%                     P4y = trace(rho_middle1*sigma_y);
%                     P4z = trace(rho_middle1*sigma_z);
%                     P4 = sqrt(P4x^2 + P4y^2 + P4z^2);
% 
% 
%                   A4 = [middle_spinn;P4x;P4y;P4z;P4];
%                   fprintf(fileID1_4,'%10.18f  \t %10.18f \t %10.18f \t %10.18f \t %10.18f\n',A4);


end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
end
 
         end
 
	% .............................................................................................................................................................
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Calculating rho for last spin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if loc_spin_last_state == 1

	fprintf('Calculating spin density for last spin of the system\n')


	% ..................................................................................


	% 1st element of the matrix |up_last><up_last| .....................................

		rho_last = sparse(2,2);

		for l1  = 1:2
			for m1 = 1:2

				rho_last(l1,m1) = sum(diag( rho_spin(l1:2:length(rho_spin),m1:2:length(rho_spin)) ));

			end
		end
	full(rho_last)

	% ..................................................................................
	% Taking observation of last spin dynamics ...........................................................		
                      A11 = [N_loc_spin; t/(0.6591);trace(rho_last*sigma_x);trace(rho_last*sigma_y);trace(rho_last*sigma_z); ((trace(rho_last*sigma_x))^2+(trace(rho_last*sigma_y))^2+(trace(rho_last*sigma_z))^2)^0.5;trace(rho_last*rho_last)];
                      fprintf(fileID3,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A11);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	end	% End of if N_loc_spin ~=0


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		%fprintf('Average Energy\n');
		%E_av = (Psi'*Hamilto*Psi) ;
	
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	
	%avg_trans(obs) = mean(conj(Psi((pos_loc_i-1)*2^(1+N_loc_spin) + 1:length(Psi))).*(Psi((pos_loc_i-1)*2^(1+N_loc_spin) + 1:length(Psi))));
	
		it_count;
		it_count = 0 ;
        
	end 	% Condition loop for taking observations at every it_count = 50 th time step 
		


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end	% end of time stepping 
%





video = VideoWriter('He.avi','Uncompressed AVI');
%video.FrameRate = 5;
open(video);
writeVideo(video,F);
close(video);

disp('Playing the movie in real space')





















