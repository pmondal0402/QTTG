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

set(gca,'FontSize',20,'ticklabelinterpreter','latex')
set(gca,'XMinorTick','off','YMinorTick','off')

% set(gca,'LineWidth',0.5,'TickLength',[0.020 0.020]);

