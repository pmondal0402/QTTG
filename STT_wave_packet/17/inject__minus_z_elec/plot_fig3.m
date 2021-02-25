M = dlmread('electron_data.txt');
t = M(1:280,1);
Px_elec =M(1:280,2);
Py_elec = M(1:280,3);
Pz_elec = M(1:280,4);
P_elec = M(1:280,5);



M1 = dlmread('Purity_1st_loc.txt');
Px_1 = M1(1:280,3);
Py_1 = M1(1:280,4);
Pz_1 = M1(1:280,5);
P_1 = M1(1:280,6);

M4 = dlmread('Purity_last_loc.txt');
Px_4 = M4(1:280,3);
Py_4 = M4(1:280,4);
Pz_4 = M4(1:280,5);
P_4 = M4(1:280,6);

M_p = dlmread('tx_Purity_elec_sub_space.txt');
purity_elec = M_p(1:280,2);

M_loc = dlmread('tx_Purity_loc_sub_space.txt');
purity_loc = M_loc(1:280,2);
x_pos =  50;
xmax = 400;
a = subplot(2,2,1)


plot(t,Px_elec,'r:','LineWidth',1.5);
hold on



%{
plot(t(1:10:end),Py_elec(1:10:end),'ro','MarkerSize',6,'MarkerFaceColor','none');
plot(t(1:1:end),Pz_elec(1:1:end),'rs','MarkerSize',6,'MarkerFaceColor','none');
%}

plot(t(1:1:end),Py_elec(1:1:end),'r-.','LineWidth',1.5);
plot(t,Pz_elec,'r--','LineWidth',1.5);


plot(t,P_elec,'r-','LineWidth',1.5);
ylabel('Polarization','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')
ylim([-0.3,1.1]);
xlim([0,xmax]);
lgnda = legend('P_x^e','P_y^e','P_z^e','P^e','Interpreter','LaTex')
set(lgnda,'color','none','FontName','Times New Roman','FontSize',10)
legend('boxoff');
pbaspect([1.3 1 1]);
text(-x_pos,1.05,'(a)','FontSize',20,'Interpreter','LaTex','FontName','Times New Roman')
hold off

b = subplot(2,2,2)
plot(t,purity_elec,'r','LineWidth',1.5);
hold on
plot(t(1:20:end),purity_loc(1:20:end),'s','MarkerSize',4,'MarkerFaceColor','blue')
ylabel('Purity','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman');
lgndb = legend('Electron Subsystem','Local Spin Subsystem','Interpreter','LaTex')
set(lgndb,'color','none','FontName','Times New Roman','FontSize',10)
legend('boxoff');
xlim([0,xmax]);
ylim([-0.3,1.1]);


pbaspect([1.3 1 1]);
text(-x_pos,1.05,'(b)','FontSize',20,'Interpreter','LaTex','FontName','Times New Roman')
hold off

c = subplot(2,2,3)

plot(t,Px_1,'b:','LineWidth',1.5);
hold on

plot(t(1:1:end),Py_1(1:1:end),'b-.','LineWidth',1.5);
plot(t(1:1:end),Pz_1(1:1:end),'b--','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','none');


plot(t,P_1,'b-','LineWidth',1.5);
ylabel('Polarization','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')

ylim([-0.3,1.1]);
xlim([0,xmax]);


pbaspect([1.3 1 1]);
lgndc = legend('P_x^1','P_y^1','P_z^1','P^1','Interpreter','LaTex')
set(lgndc,'color','none','FontName','Times New Roman','FontSize',10)

legend('boxoff')
text(-x_pos,1.05,'(c)','FontSize',20,'Interpreter','LaTex','FontName','Times New Roman')
hold off
xlabel('Time (fs)','Interpreter','LaTex','FontSize',20)


d = subplot(2,2,4)

plot(t,Px_4,'b:','LineWidth',1.5);
hold on


plot(t(1:1:end),Py_4(1:1:end),'b-.','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','none');
plot(t(1:1:end),Pz_4(1:1:end),'b--','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','none');



%{
plot(t(1:50:end),Py_4(1:50:end),'bo','MarkerSize',6,'MarkerFaceColor','none');
plot(t,Pz_4,'b-.','LineWidth',1.5);
%}
plot(t,P_4,'b-','LineWidth',1.5);



ylabel('Polarization','Interpreter','LaTex','FontSize',20,'FontName','Times New Roman')
ylim([-0.3,1.1]);
xlim([0,xmax]);

lgndd = legend('P_x^{5}','P_y^{5}','P_z^{5}','P^{5}','Interpreter','LaTex')
set(lgndd,'color','none','FontName','Times New Roman','FontSize',10)
legend('boxoff')
pbaspect([1.3 1 1]);
text(-x_pos,1.05,'(d)','FontSize',20,'Interpreter','LaTex','FontName','Times New Roman')
hold off

set(a,'LineWidth',2,'YMinorTick','on','XMinorTick','on')
set(b,'LineWidth',2,'YMinorTick','on','XMinorTick','on')
set(c,'LineWidth',2,'YMinorTick','on','XMinorTick','on')
set(d,'LineWidth',2,'YMinorTick','on','XMinorTick','on')
posa = get(a,'Position')
posb = get(b,'Position')
posc = get(c,'Position')
posd = get(d,'Position')

posb(1) = 0.44;
set(b,'Position',posb);
posd(1) = 0.44;
set(d,'Position',posd);
