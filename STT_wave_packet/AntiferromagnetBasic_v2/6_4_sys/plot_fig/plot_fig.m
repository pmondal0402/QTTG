clc; clear all; close all; 
% Spin w/o impurity
Spin = dlmread('../12Spin.txt') ;
Ntot = length(Spin) ;  % Total number of spin
% Spin with impurity
Spin_Imp = dlmread('../12SpinImp1.txt') ; 


subplot(2,1,1)
plot(Spin(:, 3), 'r-.s', 'LineWidth', 2, 'MarkerFaceColor', 'r') ; 
label(0, 1, 0, '$\mathrm{Site}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','on','YMinorTick','on')
ylim([-0.01 0.01])
xlim([1 Ntot])

subplot(2,1,2)
plot(Spin_Imp(:, 3), 'r-.s', 'LineWidth', 2, 'MarkerFaceColor', 'r') ; 
label(1, 1, 0, '$\mathrm{Site~Index}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','on','YMinorTick','on')
xlim([1 Ntot])
% Panel no. 
add_txt(11.5,0.155 , '$\mathrm{(b)}$', 'k')
add_txt(11.5,0.85 , '$\mathrm{(a)}$', 'k')


fig_layout(1)
saveas(gcf, '12spinplot.pdf')
