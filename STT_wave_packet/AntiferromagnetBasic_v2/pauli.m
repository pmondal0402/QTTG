function [reslt] = pauli(str)
% Returns pauli matrices

if str == 0
% Identity
 reslt = eye(2);
end

if str == 1
  % sigmax
  reslt = [0 1;1 0];
end

if str == 2
  % sigmay
   reslt =  [0 -1i;1i 0];
end

if str == 3
  % sigmaz
   reslt = [1 0;0 -1];
end
