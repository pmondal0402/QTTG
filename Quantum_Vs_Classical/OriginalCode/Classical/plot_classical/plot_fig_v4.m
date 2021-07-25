clc ; clear all; close all; 
% Classicl data
data_CL = dlmread('../res/debug2.txt') ;
% QM data
data = dlmread('../../res/spinQM.txt') ; 
% Electron entanglement data
path1 =...
'~/QM_transport/ClassiclaVsQM/results/UnderstandingUtkbCode/3elec1Spin/job7/06/TryFourSpins/WithSingleElec/res';  
S_en_elec = dlmread(strcat(path1, '/entropy_elec.txt')) ; 

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

ax1 = subplot(2,2,1)
% Plot S1 
% figure
plot(t_QM, data(:, 5+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st:end), data_CL(1:st:end, 5+val)/sval, 'ko')
label(0,1,0,'$\mathrm{Time~(fs)}$','$\mathrm{\langle\hat{S}_2^x\rangle~(S\hbar)}$')
yticks([-0.04:0.02:0.01])
xlim([0 xmax1*0.01])
ylim([-0.04 0.02])
set_axis_pos(ax1, 0.128, 0, 0, 0)

ax2 = subplot(2,2,2)
plot(t_QM, data(:, 6+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st:end), data_CL(1:st:end, 6+val)/sval, 'ko')
label(0,1,0,'$\mathrm{Time~(fs)}$','$\mathrm{\langle \hat{S}_2^y\rangle~(S\hbar)}$')
yticks([-0.03:0.02:0.03])
xlim([0 xmax1*0.01])
% pbaspect([1 1 1])
ylim([-0.01 0.019])
set_axis_pos(ax2, 0.60, 0, 0, 0)

ax3 = subplot(2,2,3)
plot(t_QM, data(:, 7+val), 'k-', 'LineWidth', 2)
hold on
plot(t_CL(1:st3:end), data_CL(1:st3:end, 7+val)/sval, 'ko')
hold off
label(1,1,0,'$\mathrm{Time~(fs)}$','$\mathrm{\langle \hat{S}_2^z\rangle~(S\hbar)}$')
yticks([-1.0:0.02:0.96])
xlim([0 xmax1*0.01])
% pbaspect([1 1 1])
ylim([-1 -0.96])
set_axis_pos(ax3, 0.128, 0.135, 0, 0.35)

ax4 = subplot(2,2,4)
% Get entropy data for S = 1 hbar ; 1 electron 
S1_en = dlmread('../../res/S_enHalf.txt') ; 
S1_en_thre_elec = dlmread('../../../WithThreeElec/res/S_enHalf_new.txt') ; 
yyaxis left
plot(t_QM, S1_en, '-', 'LineWidth', 2, 'color', rgb('orange'))
hold on
plot(t_QM(1:xmax1), S1_en_thre_elec(1:xmax1),...
             '-', 'LineWidth', 2, 'color', rgb('green'))
ylim([-0.16 0.16])
yyaxis right
plot(t_QM(1:xmax1), S_en_elec(1:xmax1, 2),...
             '-', 'LineWidth', 2, 'color', rgb('cyan'))
ylim([0. 0.1])
yyaxis left
label_2(1,1,0,'$\mathrm{Time~(fs)}$','Entanglement Entropy $\mathcal{S}^\mathrm{lspins}_\frac{N}{2}$','')
xlim([0 xmax1*0.01])
ax = gca ;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set_axis_pos(ax4, 0.60, 0.135, 0, 0.35)
% yticks([0:0.05:0.2])

fig_layout(1)

% add labels

xd = 0; yd = 0.18 ; 
% xa = xd+0.5; ya = 0.58 ; 
xc = xd - 112 ; yc = yd ; 
xb = xd ; yb = 0.584 ; 
xa = xc; ya = yb ; 

add_txt(xd, yd, '$\mathrm{(h)}$', 'k')
add_txt(xc, yc, '$\mathrm{(g)}$', 'k')
add_txt(xa, ya, '$\mathrm{(e)}$', 'k')
add_txt(xb, yb, '$\mathrm{(f)}$', 'k')

add_txt(xd + 20, 0.0, '$\mathrm{N_e} = 1$', rgb('orange'))
add_txt(xd + 45, 0.0, '$\mathrm{N_e} = 3$', rgb('green'))

saveas(gcf, 'AFM_Jsd_timefig1_v5.pdf')
% Note : In paper the name is changed to Fig4.pdf
% TODO : Align ylabels 
