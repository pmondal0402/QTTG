function [Spin2D] = vec_to_matrix(Spin2D_raw, Nx, Ny)
% Converts row indices of vector to 2D matrix
Spin2D = zeros(Ny, Nx) ;
pos = 0 ; 
for ii = 1:Ny
   Spin2D(ii, 1:Nx) = Spin2D_raw(1 + pos:Nx + pos) ;
   pos = pos + Nx ; 
end

 

