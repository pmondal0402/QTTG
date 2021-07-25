% NOTE : Original code has been modified to compute Hamiltonian for 2D square
% lattice

function [H, S1x] = F_aadag_ham(Jh,Nspins,n,s, Bz, Lx, Ly)
  % Call a-operator
  a = F_Aop(n);
  
  % Anisotropy 
  ani =  0*10^(-4);
  
  % Construct the ladder operators
  dim = n^Nspins;
  af = zeros(dim,dim,Nspins);
  adagf = zeros(dim,dim,Nspins);
  
  for c=1:Nspins
      nl = (c-1);
      nr = Nspins - c;
      diml = n^nl;
      dimr = n^nr;
  
      af(:,:,c) = kron(kron(eye(diml),a),eye(dimr));
      adagf(:,:,c) = kron(kron(eye(diml),a'),eye(dimr));
  end
  
  % Construct the Hamiltonian
  % TODO 

  H = zeros(dim,dim);
  for c=1:Nspins-1
  
    Spi = sqrt(2*s)*sqrtm(speye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s))*af(:,:,c);
    Smi = sqrt(2*s)*adagf(:,:,c)*sqrtm(eye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s));
    Szi = (s*speye(dim) - adagf(:,:,c)*af(:,:,c));
    
    Spii = sqrt(2*s)*sqrtm(speye(dim) - adagf(:,:,c+1)*af(:,:,c+1)/(2*s))*af(:,:,c+1);
    Smii = sqrt(2*s)*adagf(:,:,c+1)*sqrtm(speye(dim) - adagf(:,:,c+1)*af(:,:,c+1)/(2*s));
    Szii = (s*speye(dim) - adagf(:,:,c+1)*af(:,:,c+1));
   
    % Add Horizontal bonds 
    if c == Lx
       continue
    end

    H = H - Jh*0.5*(Spi*Smii + Smi*Spii) - (Jh+ani)*Szi*Szii;
     
  end
  
  % TODO : Add generalized code instead of bruitforce like below 
  % Add vertical interaction

  if Ly > 1
   for c = 1:Lx
    c ;  
    Spi = sqrt(2*s)*sqrtm(speye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s))*af(:,:,c);
    Smi = sqrt(2*s)*adagf(:,:,c)*sqrtm(speye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s));
    Szi = (s*speye(dim) - adagf(:,:,c)*af(:,:,c));
    
    Spii = sqrt(2*s)*sqrtm(speye(dim) - adagf(:,:,c+Ly)*af(:,:,c+Ly)/(2*s))*af(:,:,c+Ly);
    Smii = sqrt(2*s)*adagf(:,:,c+Ly)*sqrtm(speye(dim) - adagf(:,:,c+Ly)*af(:,:,c+Ly)/(2*s));
    Szii = (s*speye(dim) - adagf(:,:,c+Ly)*af(:,:,c+Ly));

    H = H - Jh*0.5*(Spi*Smii + Smi*Spii) - (Jh+ani)*Szi*Szii;
   end
  end
  
  % TODO : Add a generalized code for 2D square lattice for B-external
  % Manually setting external B-field for first chain i.e. Ly = 1 
  for c = 1:Lx% Nspins 
     fac = (-1)^(c+1) ;   
     H = H + fac*Bz* (s*speye(dim) - adagf(:,:,c)*af(:,:,c));
  end


  % Manually setting external B-field for first chain i.e. Ly = 2 
  for c = Lx+1:Nspins 
     fac = (-1)^c ;   
     H = H + fac*Bz* (s*eye(dim) - adagf(:,:,c)*af(:,:,c));
  end
  
  c = 1 ; 
  Spi = sqrt(2*s)*sqrtm(eye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s))*af(:,:,c);
  Smi = sqrt(2*s)*adagf(:,:,c)*sqrtm(eye(dim) - adagf(:,:,c)*af(:,:,c)/(2*s));
  S1x = 0.5*(Spi + Smi) ; 
end


