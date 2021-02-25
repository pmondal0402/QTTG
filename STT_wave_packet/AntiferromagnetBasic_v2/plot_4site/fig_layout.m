function [] = fig_layout(layout)

if layout==1
   % Fig layout
   set(gcf, 'Units','Normalized','OuterPosition',[0.25 0, 0.62, 1]) 
   set(gcf,'papersize',[12 12]);
end
