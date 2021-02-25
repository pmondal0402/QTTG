function [spinor_state] = spinor(Px,Py,Pz,P)

theta =acos( Pz/P);

if Px == 0 && Py == 0
	phi = 0;
else
  if Px == 0
     phi = asin(Py/(P*sin(theta)));
  else
     phi = acos(Px/(P*sin(theta)));
  end

end
%phi = atan(Py/Px);
spinor_state = (cos(theta/2)*[1 0]' + exp(1i*phi)*sin(theta/2)*[0 1]');
%spinor_state = spinor_state/spinor_state'*spinor_state;

return 

N = P;
nx = Px;
ny = Py;
nz = Pz;

sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

Sn = nx*sigma_x + ny*sigma_y + nz*sigma_z;
[V,E] = eig(Sn);

%spinor_state = V(:,2);

