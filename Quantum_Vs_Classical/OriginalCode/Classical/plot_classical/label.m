function [] = label(xlab, ylab, zlab, strx, stry, strz)

if xlab == 1
xlabel(strx,'Interpreter','LaTex');
end

if ylab == 1
ylabel(stry,'Interpreter','LaTex');
end

if zlab == 1
zlabel(strz,'Interpreter','LaTex');
end

set(gca,'FontSize',26,'ticklabelinterpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')

set(gca,'LineWidth',1,'TickLength',[0.020 0.020]);
xticks([0:20:100])
