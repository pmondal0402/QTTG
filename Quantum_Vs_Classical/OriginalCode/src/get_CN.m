function[Psi] = get_CN(Hamilto, Psi, dt, dirsolver, itersolver)

Id = speye(size(Hamilto));
Hp = Id + 0.5*1i*Hamilto*dt;
Hm = Id - 0.5*1i*Hamilto*dt;

% Using as solver the direct solver
if dirsolver ==1
  Psi = Hp\(Hm*Psi) ;
end

if itersolver ==1 
  brhs = Hm*Psi;
  [Psi,ierr,relres,niter] = pcg(Hp,brhs) ;
  if ierr ~= 0
    fprintf(1,'Iterative solver report: ierr = %d \n',ierr);
    fprintf(1,'   Timestep:             %d\n',it_count);
    fprintf(1,'   Rel residuum:         %e\n',relres);
    fprintf(1,'   Number of iterations: %d\n',niter);
    pause
  end
end

