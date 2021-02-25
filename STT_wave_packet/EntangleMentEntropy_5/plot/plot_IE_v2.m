clc; clear all; close all
% Date : April 26, 2020
% Author : P. Mondal
% Plots information entropy (IE)
% Fig : IE.pdf


% The data is organized as :  IE = [t/(0.6591); S_hlf ; M_{1,l}]  
M = dlmread('../IE_.txt') ; % contains data for Ne = first four data
M11 =dlmread('../IE_N11.txt') ; % contains last data for 11 Ne

data_size = length(M(:,1)) ;
M(data_size+1, :) = M11(end,:) ;

% M(isnan(M))=0 ;
subplot(1,2,1)
plot(M(:,1), M(:,2), 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r') ;
label(1, 1, 0, '$\mathrm{Number ~of ~localized ~spins ~N}$', '$\mathrm{S(1:\frac{N}{2}+1)}$','$\mathrm{E_x,q_y)}$'); 

xticks(M(:,1))
% xticklabels({num2str(M(:,1))}) ;

subplot(1,2,2)
plot(M(:,1), M(:,3), 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r') ;
label(1, 1, 0, '$\mathrm{Number ~of ~localized ~spins ~N}$', '$\mathrm{M_{(1,N)}}$','$\mathrm{E_x,q_y)}$'); 
xticks(M(:,1))
fig_layout(1)

% Add panel no.s
x1 = 10.0;y1 = 0.635;dx1 = 0.45; 
text(x1, y1,'$\mathrm{(b)}$','FontSize',20,'Interpreter','LaTex')
text(-dx1, y1 ,'$\mathrm{(a)}$','FontSize',20,'Interpreter','LaTex')

saveas(gcf, 'IE_.pdf')
% save overall data
dlmwrite('IE_manybody.txt', M)
