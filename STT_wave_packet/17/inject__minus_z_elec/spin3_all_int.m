%% 
clc
clear all

fileID1 = fopen('tuned_data_elex_z.txt','w');
%fileID3 = fopen('Last_loc.txt','w');


%%
dirsolver = 0;
itersolver = 1;
Lx = 100;
Ly=1;

lx1 = 0;
lx2 = 3;

N_loc_spin = 3 ;
pos_loc_i = 21;

one = 0;
two = 0;
Jh= .10;%1.0;
Jsd = 1.0;%1.0;
% Bz = 0;%-10.1;%5.1;
% Bz_e = 0;%Bz;%-10.1;
% B0x = 0;%3.27;

sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

t = 1;
t_so = 0;

rashbha_x = (-t*eye(2) -1i*t_so*sigma_y);
rashbha_y = (-t*eye(2) +1i*t_so*sigma_x);

%kx = 0.1;
%deltakx = 0.3;%0.3;



Tw= 0.02;%0.02;
t0=.2;
dt = 0.001;

lambda = 0.1;%0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TB Hamiltonian
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




%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

h=   -Jh*kron(eye(2*Lx*Ly),(H_heisenberg_gen_spin(N_loc_spin,1,N_loc_spin))) ...
    + H_tot_rashba_spin...
    -Jsd*H_sd(Lx,Ly,N_loc_spin,pos_loc_i);


[V,E] = eig(full(h));
P = 0;
for g = 1:5
    P = P + V(:,g)*V(:,g)';
end


%% Building |Psi> and starting to tune to ground state
X_delta_kx = 0;
for deltakx = 0.08:0.005:0.4
    X_delta_kx = X_delta_kx + 1;
        X(X_delta_kx) = deltakx;
    Y_kx = 0;
    for kx = 0.01:0.005:0.3
    Y_kx = Y_kx + 1;
    Y(Y_kx) = kx;
    Psi_ns = zeros(Lx*Ly,1);

    %txx = zeros(tmax,1);
    index=1;
    for x=1:Lx
        for y=1:Ly
            Psi_ns(index,1) = exp(1i*kx*x - deltakx^2*(x-10)^2/4);
            index=index+1;
        end
    end

    Psi_ns = Psi_ns/sqrt(Psi_ns'*Psi_ns);
    Psi_e=kron(Psi_ns(:,1),([0 1]')) ;%+  kron(Psi_ns(:,1),[0 1]');     % x_spin polarized electron wave packet
    Psi_e = Psi_e/sqrt(Psi_e'*Psi_e);


    %sparse(Psi_e*Psi_e')
    if N_loc_spin ~= 0
        %Psi = kron( kron(Psi_e,[1 1]'), [ 1 1]');    % two local spins... Up spin
        %local electrons 
        Psi = kron(Psi_e, [0 1]');
        for i = 1: N_loc_spin-1                     % N_loc_spin no. of local spins
            Psi = kron(Psi, [0 1]');
        end
        Psi=Psi/sqrt(Psi'*Psi);
    else
        Psi=Psi_e/sqrt(Psi_e'*Psi_e);
    end

    
    
     P_eigenstate(Y_kx,X_delta_kx) = (Psi'*P*Psi);
     A = [deltakx;kx;Psi'*h*Psi;(Psi'*P*Psi)] ;
     fprintf(fileID1,'%22.16f %22.16f %22.16f %22.16f\n',A);
    
    end
end
%%
[X1,Y1] = meshgrid(X,Y);
[C,h] = contour(X1,Y1,real(P_eigenstate),'ShowText','on');
xlabel('$\delta k_x$','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')
ylabel('$k_x$','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman') 
% %c

return

a= (1+lambda^2);
Hamilto = sparse((1/a)*(h-1i*lambda*(h-(Psi'*h*Psi)*sparse(eye(size(h))))));
Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;

		%Using as solver the direct solver
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



%Psi = Hp\(Hm*Psi);
Psi = Psi/sqrt(Psi'*Psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


it_count = 0;
for i = 1:10000
    t = 0.001*i;
  
    if t == 4.0
        Jsd = 0;
    end
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
h=   -Bz*kron(eye(2*Lx*Ly),(applied_field(1))) ...
    - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(1))))...
    -Jh*kron(eye(2*Lx*Ly),(H_heisenberg_gen_spin(1,1,1))) ...
    + H_tot_rashba_spin...
    -Jsd*H_sd(Lx,Ly,N_loc_spin,pos_loc_i)...
    -Bz_e*kron(eye(Lx*Ly),kron(applied_field(1),eye(2^1)));
%h
pulse(i) = B0x*gaussian;

a= (1+lambda^2);
Hamilto = (1/a)*(h-1i*lambda*(h-(Psi'*h*Psi)*eye(size(h))));
Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;

Psi = Hp\(Hm*Psi);
Psi = Psi/sqrt(Psi'*Psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlap
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;
% 
% c1 = state1'*Psi;
% c2 = state2'*Psi;
% c3 = state3'*Psi;
% c4 = state4'*Psi;
% 
% prob1(i) = conj(c1)*c1;
% prob2(i) = conj(c2)*c2;
% prob3(i) = conj(c3)*c3;
% prob4(i) = conj(c4)*c4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st element of the matrix |up_last><up_last| .....................................
rho_spin=0;
spin=1
if two == 1
    rho_spin = V(:,spin)*V(:,spin)';
		rho_last = sparse(2,2);

		for l1  = 1:2
			for m1 = 1:2

				rho_last(l1,m1) = sum(diag( rho_spin(l1:2:length(rho_spin),m1:2:length(rho_spin)) ));

			end
		end
	full(rho_last)
        txx(i) = t;
		px(i) = trace(rho_last*sigma_x)
        py(i) = trace(rho_last*sigma_y)
        pz(i) = trace(rho_last*sigma_z)

	% ..................................................................................
	% Taking observation of last spin dynamics ...........................................................		
%                       A11 = [N_loc_spin; tx(i)/(0.6591);trace(rho_last*sigma_x);trace(rho_last*sigma_y);trace(rho_last*sigma_z); ((trace(rho_last*sigma_x))^2+(trace(rho_last*sigma_y))^2+(trace(rho_last*sigma_z))^2)^0.5;trace(rho_last*rho_last)];
%                       fprintf(fileID3,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A11);
%   
% 	% ..............


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% rho = Psi*Psi';
% txx(i) = t;
% px(i) = trace(rho*sigma_x);
% py(i) = trace(rho*sigma_y);
% pz(i) = trace(rho*sigma_z);

end
subplot(2,1,1)
plot(txx,px,'b-','LineWidth',2)
hold on
plot(txx,py,'r-','LineWidth',2)
plot(txx,pz,'k-','LineWidth',2)

%plot(tx,p,'LineWidth',2)
legend('px delk 1 k 1','py','pz')

% 
% subplot(2,1,1)
% plot(txx,prob1,'r-','LineWidth',2);
% hold on
% plot(txx,prob2,'b-','LineWidth',2);
% plot(txx,prob3,'g-','LineWidth',2);
% plot(txx,prob4,'k-','LineWidth',2);
% legend('uuu','ddd','udu','udd')
return



subplot(3,2,1)
t = 0.00;
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
spy( - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(4)))))
legend('t=0.02')


subplot(3,2,2)
t = 0.2;
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
spy( - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(4)))))
legend('t = 0.2')


subplot(3,2,3)
t = 2.0;
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
spy( - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(4)))))
legend('t = 2.0')


subplot(3,2,4)
t = 4;
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
spy( - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(4)))))
legend('t = 4.0')


subplot(3,2,5)
t = 6;
gaussian = exp(-0.5*(t-t0)^2/Tw^2);
spy( - B0x*gaussian*kron(eye(2*Lx*Ly),(kron(sigma_x,eye(4)))))
legend('t = 6.0')






[Vh,Eh]= eig(full(h));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Evolution
psi = 0;
psi = kron([1 0]',kron([1 0]',[1 1]'));
psi = psi/sqrt(psi'*psi);

for n = 1:length(Eh)
    coff(n) = Vh(:,n)'*psi;
    %coff(n) = c*conj(c);
 end
% 
% full(-h)
% Eh
% Vh
% coff
psi'*h*psi;
%%
for i = 1:1500

tx(i) = 0.1*i;
psit = 0;
for n = 1:length(Eh)
psit = psit + coff(n)*Vh(:,n)*exp(-1i*tx(i)*Eh(n,n));
end
%..........................................................
psit = psit/sqrt(psit'*psit);
rho_spin = 0;
% Considering first spin dynamics in a two spin system
rho_spin= psit*psit';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if one == 1
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

		full(rho1)
		% ...............................................................
		px1(i) = trace(rho1*sigma_x);
        py1(i) = trace(rho1*sigma_y);
        pz1(i) = trace(rho1*sigma_z);
        p(i) = sqrt(px1(i)^2 + py1(i)^2 + pz1(i)^2);


		% Taking observation of 1st spin dynamics ...........................................................
		
		A1 = [N_loc_spin/N_loc_spin; tx(i)/(0.6591);trace(rho1*sigma_x);trace(rho1*sigma_y);trace(rho1*sigma_z); ((trace(rho1*sigma_x))^2+(trace(rho1*sigma_y))^2+(trace(rho1*sigma_z))^2)^0.5;trace(rho1*rho1)];
                fprintf(fileID1,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A1);

end
%..........................................................



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st element of the matrix |up_last><up_last| .....................................
if two == 1
		rho_last = sparse(2,2);

		for l1  = 1:2
			for m1 = 1:2

				rho_last(l1,m1) = sum(diag( rho_spin(l1:2:length(rho_spin),m1:2:length(rho_spin)) ));

			end
		end
	full(rho_last)
		px2(i) = trace(rho_last*sigma_x);
        py2(i) = trace(rho_last*sigma_y);
        pz2(i) = trace(rho_last*sigma_z);

	% ..................................................................................
	% Taking observation of last spin dynamics ...........................................................		
                      A11 = [N_loc_spin; tx(i)/(0.6591);trace(rho_last*sigma_x);trace(rho_last*sigma_y);trace(rho_last*sigma_z); ((trace(rho_last*sigma_x))^2+(trace(rho_last*sigma_y))^2+(trace(rho_last*sigma_z))^2)^0.5;trace(rho_last*rho_last)];
                      fprintf(fileID3,'%4.6f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n',A11);
  
	% ..............


end

%..........................................................

end
%subplot(2,1,2)
plot(tx,px1,'b-','LineWidth',2)
hold on
plot(tx,py1,'r-','LineWidth',2)
plot(tx,pz1,'k-','LineWidth',2)

%plot(tx,p,'LineWidth',2)
legend('px','py','pz')
saveas(gcf,'heisenberg.jpg')

% for m = 1:length(px1)
%     X(m) = px1(m);
%     for n = 1:length(py1)
%         Y(n) = py1(n);
%        Z(n,m) = pz1(m); 
%     end
%     
%     
% end
