clc; clear all; close all
% Date : April 26, 2020
% Author : P. Mondal
% Plots information entropy (IE)
% Fig : IE.pdf

% The data is organized as :  IE = [t/(0.6591); S1 ; Sl ; S1_l ; S1 + Sl - S1_l]  
M = dlmread('../IE_Ns_9.txt') ;
% M(isnan(M))=0 ;
subplot(2,2,1)
plot(M(:,1), M(:,4), 'r-o', 'LineWidth', 2) ;
% hold on
% Checked if I calculated mutual information entropy correctly
% plot(M(:,1), M(:,2) + M(:,3)-M(:,4),  'k--', 'LineWidth', 2)
%plot(M(:,1), M(:,2), 'k-', 'LineWidth', 2) ;
%plot(M(:,1), M(:,3), 'r-', 'LineWidth', 2) ;
%hold off
ylim([-0.1 0.7])
xlim([0 8])
label(0, 1, 0, '$\mathrm{Time ~(fs)}$', '$\mathrm{M_{1,9}}$','$\mathrm{E_x,q_y)}$'); 
subplot(2,2,2)
plot(M(1:5:end,1), M(1:5:end,2), 'bo', 'LineWidth', 2) ;
hold on
plot(M(:,1), M(:,3), 'r-', 'LineWidth', 2) ;
hold off
xlim([0 8])
ylim([-0.1 1.1])
label(0, 1, 0, '$\mathrm{Time ~(fs)}$', '$\mathrm{S_\alpha}$','$\mathrm{E_x,q_y)}$'); 
get_leg('$\mathrm{S_1}$', '$\mathrm{S_9}$', 20)

subplot(2,2,3)
plot(M(:,1), M(:,5), 'r-o', 'LineWidth', 2) ;
xlim([0 8])
ylim([-0.1 1.6])
label(1, 1, 0, '$\mathrm{Time ~(fs)}$', '$\mathrm{S_{1,9}}$','$\mathrm{E_x,q_y)}$'); 

subplot(2,2,4)
plot(M(:,1), M(:,6), 'r-o', 'LineWidth', 2) ;
xlim([0 8])
ylim([-0.1 2.1])
label(1, 1, 0, '$\mathrm{Time ~(fs)}$', '$\mathrm{S_{1,\dots, 5}}$','$\mathrm{E_x,q_y)}$'); 

fig_layout(1)
% Add panel no.s
x1 = 7.0;y1 = 0.1;dx1 = 3.5; dy1 = 3.05;
text(x1, y1,'$\mathrm{(d)}$','FontSize',20,'Interpreter','LaTex')
text(-dx1, y1 ,'$\mathrm{(c)}$','FontSize',20,'Interpreter','LaTex')
text(x1, y1 + dy1,'$\mathrm{(b)}$','FontSize',20,'Interpreter','LaTex')
text(-dx1, y1 + dy1,'$\mathrm{(a)}$','FontSize',20,'Interpreter','LaTex')

saveas(gcf, 'IE_S9.pdf')
