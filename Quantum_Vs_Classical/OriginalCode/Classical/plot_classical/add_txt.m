function [] = add_txt(x , y, str, str2)
% --------------------------
% For 6 or higher panel figs
% --------------------------
% text(x,y,str,'FontSize',12,'Interpreter','LaTex', 'color', str2,'FontName','Times New Roman')

% ------------------------
% for 4 or less panel figs
% ------------------------
text(x,y,str,'FontSize',26,'Interpreter','LaTex', 'color', str2,'FontName','Times New Roman')


