function [curr,N_num] = get_curr_N_mag(M1,M2,M3,M4,M5)

t = M1(:,2);
Sz1 = M1(:,5);
Sz2 = M2(:,4);
Sz3 = M3(:,4);
Sz4 = M4(:,4);
Sz5 = M5(:,5);

S = 5;

% -------------
% MAGNON NUMBER
% -------------
N_num = 0;
m_num = length(t);	% NUMEBR OF MEASUREMENTS
for m = 1:m_num
N_num(m,1) = S - (Sz1(m) + Sz2(m) + Sz3(m) + Sz4(m) + Sz5(m));
end

% -------
% CURRENT
% -------

% CONVERT WITH CONVERSION FACTOR LATER

for m = 1:m_num
curr(m,1) = (N_num(m,1))/(t(m));
end


