clc ; clear all; close all; 
% Classicl data
data_CL = dlmread('../res/debug2.txt') ;
% QM data
data = dlmread('../../res/spinQM.txt') ; 
% TODO change label 
t_CL = data_CL(:, 1) ; 
st = 300 ; 
st3 = 800 ; 
sval = 1 ; % hbar
t_QM = data(:, 1) ; 
val = -3 ; 
xmax1 = 8000 ; 

% 3 electron case 
data_QM_three_elec = dlmread('../../../WithThreeElec/res/spinQM.txt') ; 
data_CL_three_elec = dlmread('../../../WithThreeElec/Classical/res/debug2.txt') ;

% Plot S2 dynamics : x,y,z
val = 0 ; 

ax1 = subplot(2,3,1)
% Plot S1 
% figure
plot(t_QM, data(:, 5+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st:end), data_CL(1:st:end, 5+val)/sval, 'ko')
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^x\rangle$')
yticks([-0.04:0.02:0.01])
xlim([0 xmax1*0.01])
ylim([-0.04 0.005])
% pbaspect([1 1 1])
% set_axis_pos(ax1, 0., 0, 0, 0)

ax2 = subplot(2,3,2)
plot(t_QM, data(:, 6+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st:end), data_CL(1:st:end, 6+val)/sval, 'ko')
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^y\rangle$')
yticks([-0.03:0.02:0.03])
xlim([0 xmax1*0.01])
% pbaspect([1 1 1])
ylim([-0.01 0.019])
% set_axis_pos(ax2, 0., 0, 0, 0)

ax3 = subplot(2,3,3)
plot(t_QM, data(:, 7+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st3:end), data_CL(1:st3:end, 7+val)/sval, 'ko')
hold off
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^z\rangle$')
yticks([-1.0:0.02:0.96])
xlim([0 xmax1*0.01])
% pbaspect([1 1 1])
ylim([-1 -0.96])
% set_axis_pos(ax3, 0., 0, 0, 0)

ax4 = subplot(2, 1,2)
% Get entropy data for S = 1 hbar ; 1 electron 
S1_en = dlmread('../../res/S_enHalf.txt') ; 
S1_en_thre_elec = dlmread('../../../WithThreeElec/res/S_enHalf_new.txt') ; 
plot(t_QM, S1_en, '-', 'LineWidth', 2, 'color', rgb('orange'))
hold on
plot(t_QM(1:xmax1), S1_en_thre_elec(1:xmax1),...
             '-', 'LineWidth', 2, 'color', rgb('green'))
hold off
label_2(1,1,0,'$\mathrm{Time~(fs)}$','Entanglement Entropy $\mathcal{S}^\mathrm{lspins}_\frac{N}{2}$','')
xlim([0 xmax1*0.01])
set_axis_pos(ax4, 0, 0.15, 0, 0)


fig_layout(1)

% add labels

xd = 0; yd = 0.159 ; 
xa = xd+0.5; ya = 0.352 ; 
xb = xa + 29 ; yb = ya ; 
xc = xb + 29 ; yc = yb ; 

add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xc, yc, '$\mathrm{(c)}$', 'k')

add_txt(xd + 20, 0.03, '$\mathrm{N_e} = 1$', rgb('orange'))

add_txt(xd + 20, 0.017, '$\mathrm{N_e} = 3$', rgb('green'))

%{
xh = 0; yh = 0.159 ; 
xg = xh-125; yg = yh ; 
xf = xh - 247; yf = yh ; 
xe = xh - 372; ye = yh ; 
xd = xh; yd = 0.61 ; 
xc = xg; yc = yd ;
xb = xf; yb = yd ;
xa = xh; ya = yd ;


add_txt(xh, yh, '$\mathrm{(h)}$', 'k')
add_txt(xg, yg, '$\mathrm{(g)}$', 'k')
add_txt(xf, yf, '$\mathrm{(f)}$', 'k')
add_txt(xe, ye, '$\mathrm{(e)}$', 'k')
add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xc, yc, '$\mathrm{(c)}$', 'k')
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')
%}

% saveas(gcf, 'AFM_Jsd_timefig1_v2.pdf')
% TODO : Align ylabels 
