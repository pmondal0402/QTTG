function[avg_sp] = measr_spin(psi, J1x, J1y, J1z,...
                                   J2x, J2y, J2z,...
                                   J3x, J3y, J3z,...
                                   J4x, J4y, J4z) 
% Measure spin density for 2x2 spins 
% TODO : Generalize the code : Currently, it does things manually
avg_sp(1) = psi'*J1x*psi ;
avg_sp(2) = psi'*J1y*psi ;
avg_sp(3) = psi'*J1z*psi ;
 
avg_sp(4) = psi'*J2x*psi ;
avg_sp(5) = psi'*J2y*psi ;
avg_sp(6) = psi'*J2z*psi ;

avg_sp(7) = psi'*J3x*psi ;
avg_sp(8) = psi'*J3y*psi ;
avg_sp(9) = psi'*J3z*psi ;

avg_sp(10) = psi'*J4x*psi ;
avg_sp(11) = psi'*J4y*psi ;
avg_sp(12) = psi'*J4z*psi ;
