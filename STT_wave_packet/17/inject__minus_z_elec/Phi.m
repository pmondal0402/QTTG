function[PHI_angle] = Phi(Px,Py,Pz,P)

theta =acos( Pz/P) ;
if Px == 0 && Py == 0 
        phi_1 = 0;
else  
	norm =sqrt( Px^2 + Py^2)
	Px = Px/norm;
	Py = Py/norm
        phi_1 = atan2(Py,Px); 
end

PHI_angle = phi_1;

