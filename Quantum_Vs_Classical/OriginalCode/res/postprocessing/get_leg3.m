function [] = get_leg(str1, str2, str3, fontsz)
% Author P. Mondal
% Date : Mar 10, 2020
% Returns legend of the str1, str2
% Input : str1, str2, fonytsz

lgndb = legend(str1,str2,str3, 'Orientation','horizontal','Location','north');
set(lgndb,'Interpreter','LaTex')
set(lgndb,'color','none','FontName','Times New Roman','FontSize',fontsz)
legend('boxoff')
