% ------------------------------------------------------------------------
% Compared numerical result with analytical formula provided by branislav
% Fig produced : S_1NM.pdf
% Data produced : 
% ------------------------------------------------------------------------
 
clc; clear all; close all;
% num_s = [3;5;7;9;11] ;
num_s = 3:2:600 ;

for ii = 1: length(num_s)
 N = num_s(ii) ; 
 fac1 =  N/(2*N-1) ;
 fac2 = (N-1)/(2*N-1) ;
 S_1N(ii) = 0.0 + fac1*log((2*N-1)/N)/log(2)...
                + fac2*log((4*N-2)/(N-1))/log(2) ;
end

data(:,1) = num_s ; data(:,2) = S_1N ;
dlmwrite('IE_manybody_v2.txt', data) ;

% plot(data(:,1), data(:,2), 'b-o', 'LineWidth', 2, 'MarkerSize', 10)
plot(data(:,1), 2-data(:,2), 'b-o', 'LineWidth', 2, 'MarkerSize', 10)
label(1, 1, 0, '$\mathrm{Number ~of ~localized ~spins ~N_{FM}}$', '$\mathrm{M_{(1,N_{FM})}}$','$\mathrm{E_x,q_y)}$');
fig_layout(1)
saveas(gcf, 'M_1NFM_analyt.pdf')

return

% Checked get_entropy working fine with 3e-3s rho example
% We have analytical result.
test_rho = zeros(4,4);
test_rho(1,1) = 4/20 ;
test_rho(2,2) = 6/20 ;
test_rho(3,3) = 6/20 ;
test_rho(4,4) = 4/20 ;
test_rho(2,3) = 6/20 ;
test_rho(3,2) = 6/20 ;
get_entropy(test_rho) ;

figure(2)
plot(num_s, S_1N, 'b-o', 'LineWidth', 2, 'MarkerSize', 12)
hold on
S_1_l = [1.3710 ;   1.4355 ;   1.4573 ;   1.4681 ;   1.4746] ;
plot(num_s, S_1_l, 'r-*', 'LineWidth', 2, 'MarkerSize', 10)
hold off
xticks(num_s)
label(1, 1, 0, '$\mathrm{Number ~of ~localized ~spins ~N_{FM}}$', '$\mathrm{S_{(1,N_{FM})}}$','$\mathrm{E_x,q_y)}$');
get_leg('$\mathrm{Analytical}$', '$\mathrm{Numerical}$', 20)
fig_layout(1)
saveas(gcf, 'S_1NFM.pdf')
