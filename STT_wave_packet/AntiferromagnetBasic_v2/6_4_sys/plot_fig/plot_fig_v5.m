% Plotting for sys 6 x 4 
clc; clear all; close all; 
% Spin w/o impurity
Spin = dlmread('../12Spin.txt') ;
Ntot = length(Spin) ;  % Total number of spin
% Spin with impurity
Spin_Imp = dlmread('../12SpinImp1.txt') ; 
Sites = 1:length(Spin) ; 
col = 60 ; % number of colors for colorplot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing reverse axis
% a = (1:10)';
% b = rand(10, 1);
% c = rand(10, 1);
% figure
% h1 = axes
% bar(a, c)
% set(h1, 'Ydir', 'reverse')
% set(h1, 'YAxisLocation', 'Right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2D data
Nx = 4 ; Ny = 6 ; 
Spin2D_raw = dlmread('../test_priya_keldysh.txt') ;
Spin2D_Imp_raw = dlmread('../test_priya_impu_kohn.txt') ; 

Spin2D = vec_to_matrix(Spin2D_raw, Nx, Ny) ;

Spin2D_Imp = vec_to_matrix(Spin2D_Imp_raw, Nx, Ny) ; 

% Plot dashed line
% x(1:2) = 1.5 ; 
y(1) = -10 ; 
y(2) = 10 ;  


% Plot fig 
subplot(2,2,1)
if sum(Spin(:, 3)) < 10^(-10)
   % Since all spins at each site for pure case is 10^(-16), I am setting them
   % manually to zero since colorplot is giving some error otherwise. 
   Spin(:, 3) = 0 ; 
end
% caxis([-0.1 0.1])
% imagesc(Spin(:,3)')
imshow(Spin(:, 3)')
colormap(jet(col)); 
hold on
for num_dash = 1:length(Spin)-1 % We do not need a dash border for last site
x(1:2) = num_dash + 0.5 ; 
plot(x, y, 'k-', 'LineWidth', 0.1)
end

hold off
% colorbar
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 0, 0, '$\mathrm{Site~Index}$','$\mathrm{Site~Index}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca, 'YTick', [])
set(gca, 'XTick', [1:length(Spin)])
axis on
% h1 = colorbar
% set(h1,  'FontName', 'times new roman', 'FontSize', 20, 'Location', ...
%                     'northoutside')
% set(h1, 'Position', [0.1297    0.8838    0.3347    0.0094])


subplot(2,2,2)
% imagesc(Spin_Imp(:,3)')
imshow(Spin_Imp(:, 3)')
colormap(jet(col)); 
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 0, 0, '$\mathrm{Site~Index}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca, 'YTick', [])
set(gca, 'XTick', [1:length(Spin)])
axis on
h2 = colorbar
set(h2,  'FontName', 'times new roman', 'FontSize', 20, 'Location', ...
                    'northoutside')
h2.Position
set(h2, 'Position', [0.3600    0.8639    0.3357    0.0127])

subplot(2,2,3)

if sum(sum(Spin2D)) < 10^(-10)
   % Since all spins at each site for pure case is 10^(-16), I am setting them
   % manually to zero since colorplot is giving some error otherwise. 
   Spin2D(:, :) = 0 ; 
end
% caxis([-0.1 0.1])
imshow(Spin2D)
% imagesc(Spin(:,3)')
colormap(jet(col)); 
hold on
for num_dash = 1:Ny-1 % We do not need a dash border for last site
x(1:2) = num_dash + 0.5 ; 
plot(x, y, 'k-', 'LineWidth', 0.1)
end

x(1) = -10 ; 
x(2) = 10 ;  


for num_dash = 1:Ny-1 % We do not need a dash border for last site
y(1:2) = num_dash + 0.5 ; 
plot(x, y, 'k-', 'LineWidth', 0.1)
end

hold off
% colorbar
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 1, 0, '$\mathrm{Site~Index}$','$\mathrm{Site~Index}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')

set(gca, 'XTick', [1:Nx])
axis on
set(gca,'YDir','normal')

subplot(2,2,4)
imshow(Spin2D_Imp)
colormap(jet(col)); 
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 0, 0, '$\mathrm{Site~Index}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')

set(gca,'YDir','normal')
set(gca, 'XTick', [1:Nx])
axis on

fig_layout(1)

% Add panel number
dy = 0.24 ; 
xd = -1;   yd = 7 ;   
xc = -10.5 ; yc = yd ;
xb = xd ;    yb = yd + 6.1 ;
xa = xc ;    ya = yb ;    
add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xc, yc, '$\mathrm{(c)}$', 'k')
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')

saveas(gcf, '6_4spinplot_v5.pdf')



