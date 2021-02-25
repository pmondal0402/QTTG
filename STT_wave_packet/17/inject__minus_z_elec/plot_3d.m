

close all

% ------------
% READING DATA
% ------------

M = dlmread('Purity_1st_loc.txt');
M2 = dlmread('Purity_last_loc.txt'); 





limit = length(M);
t1 = M(1:limit,2);
Px_1 = M(1:limit,3);
Py_1 = M(1:limit,4);
Pz_1 = M(1:limit,5);
P_1 = M(:,6);

tl = M2(1:limit,2);
Px_l = M2(1:limit,3);
Py_l = M2(1:limit,4);
Pz_l = M2(1:limit,5);
P_l = M2(:,6);

x = Px_1;
y = Py_1;
z = Pz_1;

x2 = Px_l;
y2 = Py_l;
z2 = Pz_l;

figure
a = subplot(1,2,1)
[x, y, z] = sphere(25);

h = surf(x, y, z);
set(h, 'FaceColor', [0.3 0.1 0.7]);
set(h, 'FaceAlpha', 0.1);
shading interp
axis equal
hold on
%theta = linspace(bla bla)
%phi = linspace(bla bla)
Px_1 = Px_1(50:200:end);
Py_1 = Py_1(50:200:end);
Pz_1 = Pz_1(50:200:end);
t1 = t1(50:200:end);
for ii = 1:length(Px_1)
	theta_angle(ii) = theta(Px_1(ii),Py_1(ii),Pz_1(ii),P_1(ii));
	phi_angle(ii) = Phi(Px_1(ii),Py_1(ii),Pz_1(ii),P_1(ii));
end
x1 = sin(theta_angle) .* cos(phi_angle);
y1 = sin(theta_angle) .* sin(phi_angle);
z1 = cos(theta_angle);

plot3(x1, y1, z1,'r','LineWidth',1.5)
grid on
xlabel('X','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman');
ylabel('Y','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman');
zlabel('Z','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')
view(120,30)
hold off

b = subplot(1,2,2)
plot(t1,Px_1,'b:','LineWidth',1.5);
hold on

plot(t1(1:1:end),Py_1(1:1:end),'k-.','LineWidth',1.5);
plot(t1(1:1:end),Pz_1(1:1:end),'r--','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','none');

lgndc = legend('P_x^1','P_y^1','P_z^1','P^1','Interpreter','LaTex')
pbaspect([1.3 1 1]);
set(lgndc,'color','none','FontName','Times New Roman','FontSize',10)
legend('boxoff')
xlabel('Time (fs)','Interpreter','LaTex','FontSize',20)

ylabel('Polarization','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')

text(-6000,1.35,'(b)','FontSize',20,'Interpreter','LaTex','FontName','Times New Roman')
set(a,'LineWidth',2,'YMinorTick','off','XMinorTick','off')
set(b,'LineWidth',2,'YMinorTick','on','XMinorTick','on')




return
subplot(2,1,2)
[x2,y2,z2] = sphere();
h2 = surf(x2,y2,z2);
set(h2, 'FaceColor', [0.1 0.1 0.1]);
set(h2, 'FaceAlpha', 0.1); 
axis equal
hold on

x2 = Px_l;
y2 = Py_l;
z2 = Pz_l;

for ii = 1:length(M2)
        theta_angle_l(ii) = theta(Px_1(ii),Py_1(ii),Pz_1(ii),P_1(ii));
        phi_angle_l(ii) = Phi(Px_1(ii),Py_1(ii),Pz_1(ii),P_1(ii)); 
end

x2 =  cos(theta_angle_l) .* cos(phi_angle_l);
y2 = cos(theta_angle_l) .* sin(phi_angle_l);
z2 = sin(theta_angle_l);
plot3(x1, y1, z1,'r','LineWidth',2) ;
view(120,30)
hold off 












return
figure(1)
h = stem3(x,y,z,'filled')
h.Color = 'm';
h.MarkerFaceColor = 'y';
view(-10,35)

% https://www.originlab.com/doc/Tutorials/Convert-Spherical-to-XYZ


grid on

xv = linspace(min(x), max(x), 20); 
yv = linspace(min(y), max(y), 20); 
[X,Y] = meshgrid(xv, yv); 
Z = griddata(x,y,z,X,Y);
figure(2) 
surfc(X, Y, Z);
grid on

return
grid on
set(gca, 'ZLim',[0 100])
shading interp


figure
plot3(Px_1,Py_1,Pz_1);

return

% -----------------
% Sample trajectory
% -----------------
t = 0:pi/50:10*pi;
st = sin(t);
ct = cos(t);

figure
plot3(st,ct,t)


