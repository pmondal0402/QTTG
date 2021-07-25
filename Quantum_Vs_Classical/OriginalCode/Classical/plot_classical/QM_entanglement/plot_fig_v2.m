clc ; clear all; close all; 
% Plot entanglement entropy of total spin density matrix and 
% electron density matrix 
% and |P| 

% Date : electron
data_elec_sp = dlmread('../../../res/elec_spin.txt') ; 
N = 4 ; % Number of local spin sites
comp_sp = 3 ; % Number of spin components --> x, y, z 
sex = sum( data_elec_sp(:, 1:comp_sp:end), 2) ; 
sey = sum( data_elec_sp(:, 2:comp_sp:end), 2) ; 
sez = sum( data_elec_sp(:, 3:comp_sp:end), 2) ; 

se_mag = sex.^2 + sey.^2 + sez.^2 ; 
se_mag = sqrt(se_mag) ; 



% data for entanglement entropy for local spins and electrons
data_sp = dlmread('../../../res/entropy_spin.txt') ; 
data_elec = dlmread('../../../res/entropy_elec.txt') ;
% Data for local spin evolution
data_spin_evol = dlmread('../../../res/spinQM.txt') ; 
S2_mag = data_spin_evol(:, 5).^2 + data_spin_evol(:, 6).^2 ...
               + data_spin_evol(:,7).^2 ; 
S2_mag = sqrt(S2_mag) ; 

S3_mag = data_spin_evol(:, 8).^2 + data_spin_evol(:, 9).^2 ...
               + data_spin_evol(:,10).^2 ; 
S3_mag = sqrt(S3_mag) ; 

dt = 0.01 ; 
tf = 100 ; 
t = 0:dt:tf ; % [fs]
dim = length(t) ; 
st = 200 ; 

% plot(t, se_mag(1:dim))


ax1 = subplot(2,2,3)
plot(t, data_elec(1:dim, 2), 'r-', 'LineWidth', 2) ;
hold on
plot(t(1:st:dim), data_sp(1:st:dim, 2), 'bo', 'MarkerSize', 4,...
                                            'MarkerFaceColor', 'b') ;
hold off
label(1,1,0,'$\mathrm{Time~(fs)}$',...
         'Entanglement Entropy $\mathcal{S}_\mathrm{Sub}$','')
get_leg('$\mathrm{Electron}$', '$\mathrm{Local~Spins}$', 26)

ax2 = subplot(2,2,4)
% yyaxis left
plot(t, 2*se_mag(1:dim), 'r-', 'LineWidth', 2) ;
hold on
% yyaxis right
plot(t, S3_mag(1:dim), 'b-', 'LineWidth', 2) ;

label(1,1,0,'$\mathrm{Time~(fs)}$',...
     '$\mathrm{Polarization~vector~|{\bf P}|~(S\hbar)}$', '')

get_leg('$\mathrm{|\langle\hat{\bf s}_e\rangle|}$',... 
              '$\mathrm{|\langle\hat{\bf S}_2\rangle|}$', 20)
ylim([0.95 1])

% Adjust axis pos

set_axis_pos(ax1, 0.1302, 0, 0, 0)
set_axis_pos(ax2, 0.60, 0, 0, 0)

xb = 0; yb = 1.003 ;
xa = -140; ya = yb ; 
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')

fig_layout(1)
% saveas(gcf, 'AFM_Entanglement.pdf')
