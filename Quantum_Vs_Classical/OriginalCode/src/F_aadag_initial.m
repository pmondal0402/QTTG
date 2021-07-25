
function psi = F_aadag_initial(theta,phi,Nspins,n,s)
psi = 1;
for a=1:Nspins
    psi1 = cos(theta(1,a)/2);
    psi2 = exp(1i*phi(1,a))*sin(theta(1,a)/2);
    state = zeros(n,1);
    state(1,1) = psi1; state(2,1) = psi2;
    psi = kron(psi,state);
end
end
