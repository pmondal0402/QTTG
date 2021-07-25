import numpy as np
from td2d.constants import UNIT, SIG_X, SIG_Y, SIG_Z

def HTB(Lx, Ly, gam):
    # Retirns Tightbinding Hamiltonian for 1D
    dim = Lx*Ly
    Ham = np.zeros((dim, dim))
    if Ly == 1:
      for i in range(Lx-1):
        Ham[i, i+1] = -gam
    Ham = Ham + np.conjugate(Ham).T
    return Ham

def Hsd(Lx, Ly, pos_sp, Mag):
   # Returns sd Hamiltonian gievn position of classical spin 
   # Lx, Ly : system geometry
   # pos_sp : positions of classical spins , list object
   # Mag : Magnetizations of entire system, tuple object

   dim = Lx*Ly
   Ham_sd = np.zeros((2*dim, 2*dim))

   # Get positions of classical spins
   for item in pos_sp:
     Ham_sd = Ham_sd + np.kron(np.dot(sitevec(Lx, Ly, item), 
                       sitevec(Lx, Ly, item).T), 
                       (SIG_X*Mag[item][0] + SIG_Y*Mag[item][1]
                       + SIG_Z*Mag[item][2]))
   
   Ham_sd = 0.5*(Ham_sd + np.conjugate(Ham_sd).T)
   return Ham_sd

def sitevec(Lx, Ly, pos):
   dim = Lx*Ly   
   site = np.zeros((dim, 1))
   site[pos, 0] = 1
   return site
   



