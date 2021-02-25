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

%% #############################################################################


% Declaring initial parameters i.e. length of chain, wavepacket parameters  etc.............

Lx = 100;
Ly = 1;
N_loc_spin = 8;     % No. of local spins in System... 
pos_loc_i = 51;%25*Ly+1;%1;%10*Ly+1;%15*Ly+1;%24*Ly + 1;%65*Ly + 1; %13*Ly + 1;%70*Ly + 1 ;      % Position where local spin starts

dt = 0.01; % Iterative solver time step
tmax = 90;
% Tight binding and rashba hopping rates .........

t = 1;

t_so = 0.;%0.1;
lx1 = 0;%11;        % where Rashbha starts
lx2 = 4;%41;        % where Rashbha ends

a = 1 ; %not sure about this
%kx = 1.5;%pi/(2*a);%2.6/a;
kx = 0.1;%pi/2;%0.0855;% 0.105;%0.0055;
%deltakx = 2./a;
deltakx = 0.1;%0.2;% 0.245;%0.0805;
Jsd = - 0.1;%-0.5;%*t;%2;%1.0;
Jheis = -0.1;%*t;%0.1;%0.1;
Jx = Jheis;
Jy = Jheis;
Jz = Jheis;
Bz = -.0;
lambda = 0.;
%Bz = -5.1;


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



%% Tight Binding Hamiltonian 

%% Define position space hamiltonians.... Upper half..TB Hamiltonian

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


H1_tot_i = kron(H1_pos_i,sparse(rashbha_x));
H1_tot_o = kron(H1_pos_o,sparse((-t*eye(2))));
H1 = H1_tot_i + H1_tot_o;


H0_tot_i = kron(H0_pos_i,sparse(rashbha_y));
H0_tot_o = kron(H0_pos_o,sparse(-t*eye(2)));
H0 = H0_tot_i + H0_tot_o;

H_tot = H0+ H1;

HH = H_tot + ctranspose(H_tot);
H_tot_rashba_spin = kron(HH,eye(2^N_loc_spin));


%% Building |Psi>

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
Psi_e=kron(Psi_ns(:,1),([1 1]')) ;%+  kron(Psi_ns(:,1),[0 1]');     % x_spin polarized electron wave packet
Psi_e = Psi_e/sqrt(Psi_e'*Psi_e);


%sparse(Psi_e*Psi_e')
if N_loc_spin ~= 0
    %Psi = kron( kron(Psi_e,[1 1]'), [ 1 1]');    % two local spins... Up spin
    %local electrons 
    Psi = kron(Psi_e, [0 1]');
    for i = 2: 4%N_loc_spin-1                     % N_loc_spin no. of local spins
        Psi = kron(Psi, [0 1]');
    end
    Psi=Psi/sqrt(Psi'*Psi);

    for i = 5: 8%N_loc_spin-1                     % N_loc_spin no. of local spins
        Psi = kron(Psi, [1 0]');
    end
    Psi=Psi/sqrt(Psi'*Psi);


else
    Psi=Psi_e/sqrt(Psi_e'*Psi_e);
end

Psi0 = Psi;
%Ly

h=  Jheis*kron(eye(2*Lx*Ly),(H_heisenberg_gen_spin(N_loc_spin,1,N_loc_spin))) ...
    + H_tot_rashba_spin...
    + Jsd*H_sd(Lx,Ly,N_loc_spin,pos_loc_i) 
    + Bz*kron(speye(2*Lx*Ly), applied_field(N_loc_spin)) ;%... 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a= (1+lambda^2);
Hamilto = sparse((1/a)*(h-1i*lambda*(sparse(h)-(Psi'*h*Psi)*(speye(size(h))))));

Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;

t = 0;
it_count = 0;
obs = 1;
fprintf(1,'Timestepping loop BGN\n');

%%
set(gca,'nextplot','replacechildren')
disp('Creating the movie')

nFrames = 200;
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);





while t <= 4*tmax 
    
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

if it_count == 100
    Psi = Psi/sqrt(Psi'*Psi);
    
   %fprintf('making movie') 
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Psi_movie = zeros(Lx,1);
    
    for j = 1:Lx
        Psi_movie(j,1) = Psi(1+(j-1)*spin_conf: j*spin_conf,1)'...
                        *Psi(1+(j-1)*spin_conf: j*spin_conf,1);
    end
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    site = (1:Lx);
    
    a = plot(site,Psi_movie,'b','LineWidth',1.5);
    
    ylabel('Probability','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')

ylim([ 0,1]);
xlim([0,Lx]);


pbaspect([1.3 1 1]);
lgndc = legend(['t =' num2str(t/0.6591)],'Interpreter','LaTex');
set(lgndc,'color','none','FontName','Times New Roman','FontSize',10)

legend('boxoff')
xlabel('Site','Interpreter','LaTex','FontSize',20)
set(a,'LineWidth',1.5)
    
    
    
    
    
    
    F(obs) = getframe(gcf);
    s = size(F(obs).cdata);

   it_count = 0;
   obs = obs + 1;
    
end
                
                
                
    
end

video = VideoWriter('He.avi','Uncompressed AVI');
%video.FrameRate = 5;
open(video);
writeVideo(video,F);
close(video);

disp('Playing the movie in real space')

%movie(F)

%movie2avi(F,'my_movie.mp4')






























return


























%% Define initial parameters and matrices
%format long
%fileID = fopen('tx_350_.txt','w');


Lx = 160;
Ly = 31;


a=1.0;
kx = .44/a;
deltakx = .1/a;

lx1 = 31;        % where Rashbha starts
lx2 = 130;        % where Rashbha ends

t = 1;
t_so = .1*t;
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

sigmax = [0 1;1 0];
sigmay = [0 -1i;1i 0];
sigmaz = [1 0;0 -1];

rashbha_x = (-t*eye(2) -1i*t_so*sigma_y);
rashbha_y = (-t*eye(2) +1i*t_so*sigma_x);

%% Define position space hamiltonians.... Upper half
% H1................................................

H1_pos = zeros(Lx*Ly,Lx*Ly);
H1_pos_i = zeros(Lx*Ly,Lx*Ly);
H1_pos_o = zeros(Lx*Ly,Lx*Ly);


for n = 1:1:((Lx-1)*Ly)
       H1_pos(n,n+Ly) = t;
      
% Fix me .................................     
    if n > (lx1-1)*Ly && n<(lx2*Ly+1) %((lx2+1)*Ly+1)
        H1_pos_i(n,n+Ly) = H1_pos(n,n+Ly);
        H1_pos_o(n,n+Ly) = 0;
    else
        H1_pos_i(n,n+Ly) = 0;
        H1_pos_o(n,n+Ly) = H1_pos(n,n+Ly);
    end    
     
       
end    
H1_pos;



% H0................................................


H0_pos = zeros(Lx*Ly,Lx*Ly);
H0_pos_i = zeros(Lx*Ly,Lx*Ly);
H0_pos_o = zeros(Lx*Ly,Lx*Ly);


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






%% Defining Rashbha region and shifting Hamiltonian to spin space


H1_tot_i = kron(H1_pos_i,sparse(rashbha_x));
H1_tot_o = kron(H1_pos_o,sparse((-t*eye(2))));
H1 = H1_tot_i + H1_tot_o;


H0_tot_i = kron(H0_pos_i,sparse(rashbha_y));
H0_tot_o = kron(H0_pos_o,sparse(-t*eye(2)));
H0 = H0_tot_i + H0_tot_o;

H_tot = H0+ H1;

H = H_tot + ctranspose(H_tot);

[eigVec, eigVal] = eig(full(H));
%[eigVec, eigVal] = eigs((H));



%% Building |Psi>

Psi_ns = zeros(Lx*Ly,1);

index=1;
for x=1:Lx
    for y=1:Ly
        Psi_ns(index,1) = sin(pi*(y)/(Ly+1)/a)*exp(1i*kx*x - deltakx^2*(x-5)^2/4);
        index=index+1;
    end
end
Psi_ns = Psi_ns/sqrt(Psi_ns'*Psi_ns);
Psi=kron(Psi_ns(:,1),[1 0]');
Psi=Psi/(Psi'*Psi);

cof = zeros(1,2*Lx*Ly);

for n=1:length(eigVal)
    cof(n) = (eigVec(:,n)'*Psi);
end

%How to calc |Psi(t)>

set(gca,'nextplot','replacechildren')
disp('Creating the movie')

%% .............................................
nFrames = 200;
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);


%..............................................


%%

for g=1:100
%clf
Psi_t = zeros(length(Psi),1);
tx = 0.1*(g-1);          % Changing time in 0.1 fs step i.e. dividing by .6591 
%tx= 2*g;
%tx = g

%tx = 1;
for n=1:length(eigVal)
    Psi_t = Psi_t + cof(n)*exp(-1i*eigVal(n,n)*tx/(0.6591))*eigVec(:,n);
end

ind = 1;
for x = 1:Lx
    for y = 1:Ly
        ind;
        Psi_t_x = [Psi_t(ind) Psi_t(ind+1)].';
        Psi_rho(x,y) = Psi_t_x'*Psi_t_x;
        ind = ind +2;
    end    
end    

surface(Psi_rho','EdgeColor','none');
axis = ([0 160 0 31]);
F(g) = getframe(gca);

%{
Psi_rho = zeros(Lx*Ly,1);
ind = 1;
for site = 1:2:2*Lx*Ly
    site;
    Psi_t_x = [Psi_t(site) Psi_t(site+1)].';
    Psi_rho(ind) = Psi_t_x'*Psi_t_x;
    ind = ind+1;
end    
%}


%plot(x,Psi_rho);

%xlabel('site');
%ylabel('rho');

s = size(F(g).cdata);
fprintf('%d %d\n', s(2), s(1))
%movie2avi(F,'my_movie.mp4')

%Psi_t_transps = Psi_t';

%% Making movie for Psi_density 

%Psi_rho =  Psi_t'*Psi_t



%{
%% Calculating density matrix in spin space i.e. rho_s which is 2x2


rho_s = zeros(2);

for x = 1:2:length(eigVal)
   rho_s(1,1) = rho_s(1,1) + Psi_t(x,1)*Psi_t_transps(1,x); 
   rho_s(1,2) = rho_s(1,2) + Psi_t(x,1)*Psi_t_transps(1,x+1);
   rho_s(2,1) = rho_s(2,1) + Psi_t(x+1,1)*Psi_t_transps(1,x);
   rho_s(2,2) = rho_s(2,2) + Psi_t(x+1,1)*Psi_t_transps(1,x+1);
    
end    

A = [tx; rho_s(1,1); real(rho_s(1,2)); imag(rho_s(1,2))];
tx
rho_s

fprintf(fileID,'%22.16f %22.16f %22.16f %22.16f\n',A);
%}

end



disp('Playing the movie in real space')

movie(F)

%save F
%movie2avi(F,'my_movie_fast.avi')
%fclose(fileID);



%% Valid only for tight binding model....Use rho matrix for calculation in general cases



%{

%% Note P_z_rho and P_z are same.... two different methods


%Making a matrix in XY space to plot
index=1;
for x=1:Lx
    for y=1:Ly
        S = [Psi_t(index) Psi_t(index+1)].';
        Spin(x,y) = S'*sigmaz*S;
        S_x(x,y) = S'*sigmax*S;
        S_y(x,y) = S'*sigmay*S;

        index=index+2;
    end
end

P_z = sum(sum(Spin));
P_x = sum(sum(S_x));
P_y = sum(sum(S_y));
%P(g) = sqrt(P_x(g)^2 + P_y(g)^2+P_z(g)^2);

%surface(Spin);
%pause(.1);
%end

%}

