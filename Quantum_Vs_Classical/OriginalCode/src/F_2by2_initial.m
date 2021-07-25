function psi = F_2by2_initial(theta,phi,Nspins)
psi = 1;
for a=1:Nspins
    psi1 = cos(theta(1,a)/2);
    psi2 = exp(1i*phi(1,a))*sin(theta(1,a)/2);
    psi = kron(psi,[psi1;psi2]);
end
end
