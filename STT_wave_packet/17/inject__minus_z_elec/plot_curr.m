close all
clc;
clear all;


% ------------------------
% READING DATA FOR ELEC -z
% ------------------------

M1 = dlmread('Purity_1st_loc.txt');
M2 = dlmread('Purity_2nd_loc.txt');
M3 = dlmread('Purity_3rd_loc.txt');
M4 = dlmread('Purity_4th_loc.txt');
M5 = dlmread('Purity_last_loc.txt');

[curr,N_num] = get_curr_N_mag(M1,M2,M3,M4,M5);
%
%subplot(2,1,1)
%plot(curr(1:400),N_num(1:400))
%title('-z elec')
%



% ------------------------
% READING DATA FOR ELEC z
% ------------------------


M11 = dlmread('Purity_1st_loc1.txt');
M22 = dlmread('Purity_2nd_loc1.txt');
M33 = dlmread('Purity_3rd_loc1.txt');
M44 = dlmread('Purity_4th_loc1.txt');
M55 = dlmread('Purity_last_loc1.txt');

[curr1,N_num1] = get_curr_N_mag(M11,M22,M33,M44,M55);


plot(curr1(1:300),N_num1(1:300),'-r','LineWidth',3);
hold on
plot(curr(1:400),N_num(1:400))

title('z elec')
ylim([-0.0 0.05])
hold off


