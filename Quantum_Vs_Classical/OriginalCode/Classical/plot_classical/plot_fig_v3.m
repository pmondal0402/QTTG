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

ax1 = subplot(2,4,1)
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

ax2 = subplot(2,4,2)
plot(t_QM, data(:, 6+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st:end), data_CL(1:st:end, 6+val)/sval, 'ko')
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^y\rangle$')
yticks([-0.03:0.02:0.03])
xlim([0 xmax1*0.01])
% pbaspect([1 1 1])
ylim([-0.01 0.019])
% set_axis_pos(ax2, 0., 0, 0, 0)

ax3 = subplot(2,4,3)
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

ax4 = subplot(2,4,4)
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
%set_axis_pos(ax4, 0, 0.15, 0, 0)


fig_layout(1)

% add labels

xd = 0; yd = 0.159 ; 
xa = xd+0.5; ya = 0.352 ; 
xc = xd - 125 ; yc = yd ; 
xb = xc -122 ; yb = yd ; 
xa = xb - 122; ya = yd ; 

add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xa, yd, '$\mathrm{(a)}$', 'k')
add_txt(xb, yd, '$\mathrm{(b)}$', 'k')
add_txt(xc, yd, '$\mathrm{(c)}$', 'k')

add_txt(xd + 20, 0.017, '$\mathrm{N_e} = 1$', rgb('orange'))
add_txt(xd + 20, 0.007, '$\mathrm{N_e} = 3$', rgb('green'))

% saveas(gcf, 'AFM_Jsd_timefig1_v3.pdf')
% TODO : Align ylabels 
