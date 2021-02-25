clc; clear all; close all; 
Spin = dlmread('../6Spin.txt') ;
Entropy = dlmread('../6Spin_entropy.txt') ; 

subplot(2,1,1)
plot(Spin(:, 1), 'r-', 'LineWidth', 2.5) ; 
hold on
plot(Spin(:, 2), 'bo', 'MarkerSize', 8) ; 
plot(Spin(:, 3), 'k-.', 'LineWidth', 1.7) ; 
label(1, 1, 0, '$\mathrm{Site}$','$\mathrm{S_i^\alpha}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','on','YMinorTick','on')
get_leg3('x','y','z',20)
ylim([-0.01 0.01])

subplot(2,1,2)
bar(Entropy,'FaceColor', 'r' )
label(1, 1, 0, '$\mathrm{Site}$','$\mathrm{Local~Entropy}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','on','YMinorTick','on')

fig_layout(1)
saveas(gcf, '6spinplot.pdf')
