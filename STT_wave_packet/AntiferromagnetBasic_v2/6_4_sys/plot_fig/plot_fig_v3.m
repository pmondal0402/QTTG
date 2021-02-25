clc; clear all; close all; 
% Spin w/o impurity
Spin = dlmread('../4Spin.txt') ;
Ntot = length(Spin) ;  % Total number of spin
% Spin with impurity
Spin_Imp = dlmread('../4SpinImp1.txt') ; 
Sites = 1:length(Spin) ; 

% 2D data
Nx = 4 ; Ny = Nx ; 
Spin2D_raw = dlmread('../../QUIMB/ex1_v2/spin2d.txt') ;
Spin2D_Imp_raw = dlmread('../../QUIMB/ex1_v2/spin2d_imp.txt') ; 

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
colormap(jet(numel(Spin(:,3)' ))); 
hold on
for num_dash = 1:Nx-1 % We do not need a dash border for last site
x(1:2) = num_dash + 0.5 ; 
plot(x, y, 'k-', 'LineWidth', 0.1)
end

hold off
% colorbar
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 1, 0, '$\mathrm{Site~Index}$','$\mathrm{Site~Index}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca, 'YTick', [1])
set(gca, 'XTick', [1:Nx])
axis on

subplot(2,2,2)
% imagesc(Spin_Imp(:,3)')
imshow(Spin_Imp(:, 3)')
colormap(jet(numel(Spin_Imp(:, 3)'))); 
caxis([-0.1 0.1])
% h2 = colorbar
% set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 0, 0, '$\mathrm{Site~Index}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca, 'YTick', [1])
set(gca, 'XTick', [1:Nx])
axis on

subplot(2,2,3)

if sum(sum(Spin2D)) < 10^(-10)
   % Since all spins at each site for pure case is 10^(-16), I am setting them
   % manually to zero since colorplot is giving some error otherwise. 
   Spin2D(:, :) = 0 ; 
end
% caxis([-0.1 0.1])
imshow(Spin2D)
% imagesc(Spin(:,3)')
colormap(jet(numel(Spin(:,3)' ))); 
hold on
for num_dash = 1:Nx-1 % We do not need a dash border for last site
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

% set(gca, 'YTick', [1:Ny])
set(gca, 'XTick', [1:Nx])
set(gca,'YtickLabel',Ny:-1:1)
axis on


ax4 = subplot(2,2,4)
imshow(Spin2D_Imp)
colormap(jet(numel(Spin_Imp(:, 3)'))); 
caxis([-0.1 0.1])
h2 = colorbar
set(h2,  'FontName', 'times new roman', 'FontSize', 20)
label(1, 0, 0, '$\mathrm{Site~Index}$','$\mathrm{S_i^z}$','$\mathrm{E_x,q_y)}$');  
set(gca,'XMinorTick','off','YMinorTick','off')
% set(gca, 'YTick', [Ny:-1:1])
set(gca, 'XTick', [1:Nx])
% set(ax4, 'YDir','reverse')
set(gca,'YtickLabel',Ny:-1:1)
axis on

fig_layout(1)

% Add panel number
dy = 0.24 ; 
xd = -0.3;   yd = 0.7-dy ;   
xc = -6.2 ; yc = yd ;
xb = xd ;    yb = -3.3-dy-0.5 ;
xa = xc ;    ya = yb ;    
add_txt(xd, yd, '$\mathrm{(d)}$', 'k')
add_txt(xc, yc, '$\mathrm{(c)}$', 'k')
add_txt(xb, yb, '$\mathrm{(b)}$', 'k')
add_txt(xa, ya, '$\mathrm{(a)}$', 'k')

saveas(gcf, '4spinplot_v3.pdf')



