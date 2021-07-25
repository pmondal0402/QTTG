% Plotting for sys 6 x 4 
clc; clear all; close all; 

% data
data_len = 501 ; % length of data to plot
proj_neel = dlmread('proj_eigstate_Neel.txt') ; 

proj_neel = proj_neel(:, 1:20) ; 
% proj_neel = proj_neel(1:data_len, :) ; 
Stotz = dlmread('SpinSz.txt') ; 
Stotz = Stotz(1:data_len, :) ; 
st = 0.001 ; 
Bval = 0.:st:10 ; 
inv_st = 100 ; 

col = 500 ; 
subplot(2,2,1)
imagesc(Stotz)
colormap(jet(col));
% caxis([0.9 1])
h1 = colorbar
set(gca,'YDir','normal')
label(1, 1, 0, '$\mathrm{Spin~Index}$','$\mathrm{External~Field~\mu_B B_{0}~(eV)}$',...
                                            '$\mathrm{E_x,q_y)}$');  
set(h1,  'FontName', 'times new roman', 'FontSize', 20)
h1.Title.String = 'S_i^z';
ticy = [1:inv_st:length(Bval)] ; 
set(gca, 'YTick', ticy) ; 
set(gca, 'YTicklabels', Bval(ticy)) ;

subplot(2,2,2)
imagesc(proj_neel)
colormap(jet(col));
set(gca,'YDir','normal')
label(1, 0, 0, '$\mathrm{Eigenstate~n}$','$\mathrm{Site~Index}$','$\mathrm{E_x,q_y)}$');  
h2 = colorbar
set(h2,  'FontName', 'times new roman', 'FontSize', 20)
h2.Title.String = 'P_n';

ticy = [1:2*inv_st:length(Bval)] ; 
set(gca, 'YTick', ticy) ; 
set(gca, 'YTicklabels', Bval(ticy)) ;
% get tick values for xaxis
ticx = [1:5:201] ; 
eig_num = 0:1:201 ; 
set(gca, 'XTick', ticx) ; 
set(gca, 'XTicklabels', eig_num(ticx)) ; 

% set(gca, 'YTick', ticy) ; 
fig_layout(1)

% Add labels
xb = 1 ; yb = 1050 ; xa = xb - 31 ;  
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, yb, '$\mathrm{(a)}$', 'k')

saveas(gcf, 'AFM_GS_B-field.pdf')

