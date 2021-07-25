clc ; clear all; close all; 
% Classicl data
data_QM = dlmread('../res/debug2.txt') ;
% QM data
data = dlmread('../../res/spinQM.txt') ; 
% TODO change label 
t_QM = data_QM(:, 1) ; 
st = 300 ; 
st3 = 800 ; 
sval = 1 ; % hbar
t = data(:, 1) ; 
val = -3 ; 


% Plot first spin dynamics
% S1x : Matching
%{
plot(t, data(:, 5+val), 'b-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 5+val)/sval, 'bo')
label(0,1,0,'','$\langle \hat{S}_1^\alpha\rangle$', '')
fig_layout(1)
%}

% S1y : Matching
%{
plot(t, data(:, 6+val), 'r-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 6+val)/sval, 'ro')
label(0,0,0,'','$\mathrm{Tr}[\rho_\mathrm{gs}|\psi_n\rangle\langle \psi_n|]$')
fig_layout(1)
%}

% S1z : QM and Cl has discrepancy
%{
plot(t, data(:, 7+val), 'k-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 7+val)/sval, 'ko')
hold off
label(0,0,0,'','$\mathrm{Tr}[\rho_\mathrm{gs}|\psi_n\rangle\langle \psi_n|]$')
fig_layout(1)
%}

% Plot second spin dynamics
% S2x : Small discrepancy
%{
plot(t, data(:, 5+val), 'b-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 5+val)/sval, 'bo')
label(0,1,0,'','$\langle \hat{S}_1^\alpha\rangle$', '')
fig_layout(1)
%}


% S2y : Discrepancy
%{
plot(t, data(:, 6+val), 'r-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 6+val)/sval, 'ro')
label(0,0,0,'','$\mathrm{Tr}[\rho_\mathrm{gs}|\psi_n\rangle\langle \psi_n|]$')
fig_layout(1)
%}


% S1z : QM and Cl has discrepancy
%{
plot(t, data(:, 7+val), 'k-')
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 7+val)/sval, 'ko')
hold off
label(0,0,0,'','$\mathrm{Tr}[\rho_\mathrm{gs}|\psi_n\rangle\langle \psi_n|]$')
fig_layout(1)
%}

%{
subplot(2,3,1)
plot(t_QM, -sval_qm/sval)
hold on
plot(t_QM, -sval_qm0)
hold off
%}

% Plot dynamics S1
subplot(2,4,1)
% Plot S1 
% figure
plot(t, data(:, 5+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 5+val)/sval, 'ko')
label(0,1,0,'','$\langle \hat{S}_1^x\rangle$', '')
yticks([0:0.1:0.3])

subplot(2,4,2)
plot(t, data(:, 6+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 6+val)/sval, 'ko')
label(0,1,0,'','$\langle \hat{S}_1^y\rangle$', '')
yticks([-0.1:0.1:0.1])

subplot(2,4,3)
plot(t, data(:, 7+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 7+val)/sval, 'ko')
hold off
label(0,1,0,'','$\langle \hat{S}_1^z\rangle$', '')
yticks([0.96:0.02:1.02])

% Plot S2 dynamics : x,y,z
val = 0 ; 

subplot(2,4,5)
% Plot S1 
% figure
plot(t, data(:, 5+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 5+val)/sval, 'ko')
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^x\rangle$')
yticks([-0.04:0.02:0.01])

subplot(2,4,6)
plot(t, data(:, 6+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st:end), data_QM(1:st:end, 6+val)/sval, 'ko')
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^y\rangle$')
yticks([-0.03:0.02:0.03])

subplot(2,4,7)
plot(t, data(:, 7+val), 'k-', 'LineWidth', 2)
hold on
plot(t_QM(1:st3:end), data_QM(1:st3:end, 7+val)/sval, 'ko')
hold off
label(1,1,0,'$\mathrm{Time~(fs)}$','$\langle \hat{S}_2^z\rangle$')
yticks([-1.0:0.02:0.96])

subplot(2, 4,4)
% Get entropy data for S = 1 hbar
S1_en = dlmread('../../res/S_enHalf.txt') ; 
plot(t, S1_en, 'k-', 'LineWidth', 2, 'color', rgb('orange'))
label(0,1,0,'$\mathrm{Time~(fs)}$','$\mathcal{S}_\mathrm{Half}$','')

subplot(2, 4,8)
% Get entropy data for S = 1 hbar
Neg_S1 = dlmread('../../res/negativity.txt') ; 
plot(t, Neg_S1(:,2), '-', 'color', rgb('red'), 'LineWidth', 2)
label(1,1,0,'$\mathrm{Time~(fs)}$','$\mathrm{Logarithmic~Negativity}$','')
yticks([0.1:0.2:0.5])
fig_layout(1)

% add labels
xh = 165; yh = 0.52 ; 
xg = -148; yg = yh ; 
xf = -458; yf = yh ; 
xe = -778; ye = yh ; 
xd = xh; yd = 1.217 ; 
xc = xg; yc = yd ;
xb = xf; yb = yd ;
xa = xe; ya = yd ;

add_txt(xh, yh, '$\mathrm{(h)}$', 'k')
add_txt(xg, yg, '$\mathrm{(g)}$', 'k')
add_txt(xf, yf, '$\mathrm{(f)}$', 'k')
add_txt(xe, ye, '$\mathrm{(e)}$', 'k')
add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xc, yc, '$\mathrm{(c)}$', 'k')
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')

% saveas(gcf, 'AFM_Jsd_timefig1_v1.pdf')
% TODO : Align ylabels 
