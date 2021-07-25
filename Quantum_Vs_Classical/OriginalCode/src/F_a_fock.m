function ret = F_a_fock(c,Nspins,n)
    
a = F_Aop(n);

    nl = (c-1);
    nr = Nspins - c;
    diml = n^nl;
    dimr = n^nr;

    ret = kron(kron(eye(diml),a),eye(dimr));

end
