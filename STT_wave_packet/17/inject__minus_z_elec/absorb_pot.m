function[H_absorb] = absorb_pot(Lx,width,r11,r22,c1)

% -------------------
%ABSORBING POTENTIAL
% -------------------
c1 = 55;        %Absorbing potential centered at x = c  
a1 = (1 - 16/c1^3);
b1 = (1 - 17/c1^3)*(1/c1^2);
r11 = 70;
r22 = 10 ;
del_k_min = 0.01;%c1/(2*(r22 - r11));
E_minim = 1;

for sitee = 1:Lx
%	xx(sitee) = c1 -5 + 0.1*sitee;
	 xx(sitee) =  c1 + 2*del_k_min*(sitee - r11); 
	 yy(sitee) = a1*xx(sitee) - b1*(xx(sitee))^3 + 4/(c1 - xx(sitee))^2 ...
			-4/(c1 + xx(sitee))^2;
end
site = (1:Lx);
figure()
plot(site,yy,'b','LineWidth',2);

figure()
plot(xx,yy,'r','LineWidth',2)

for n = 1:Lx 
 H_absorb(n,n) = E_minim*yy(n);
end

