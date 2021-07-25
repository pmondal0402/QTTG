# Implement local spin evolution with magnetic pulse along x-dir 

import numpy as np
from numpy.linalg import norm as lanorm
from scipy.linalg import expm
from td2d.constants import SIG_X, SIG_Y, SIG_Z, HBAR
import Hamiltonians as hm
import rk4density 
from numpy import linalg as LA
import llg_evol as llge


# Overall plan 
# Params :
Lx = 4
Ly = 1

# Number of electronic spins 
Ne = 1 
# Hopping for tight binding 
gam = 1 # [eV] 

# Params for sd Hamiltonian
# Number of classical spins
Nsp = Lx*Ly 
# Jsd coupling 
Jsd = +0.01 # [eV] # Change sign to -ve

# Exchange coupling for AFM
# +ve for AFM and -ve for FM

Jex = 0.01 # [eV]
# Positions of classical spins
pos_sp = [0, 1, 2, 3] # Choosing second site 
# Vals for classical spins
# Mag = [[0, 0, 1], [0, 0, 1], [0, 0, 1]]

Mag = [[0, 0, 1], [0, 0, -1], [0, 0, 1], [0, 0, -1]]

# Pulse params
t0 = 50 
Tw = 15
# B0x = -0.3
B00 = +0.01
NoB0x = 0
# Strength of external B-field
B0z = +0.06

# Which method to use for time evolution
# 1 : rk4 method, 2 : Evolve operator
method = 1 #  
# Evolve LLG
# Options : 'yes', 'no'
run_LLG = 'yes'
gfac = 2 

# rk4 params
dt = 0.01 # [fs]
tmax = 200# 20000.0  # [fs]

# Create rho in energy basis with half filled states
# i.e. |1up 2up 3up> for electrons
rhoE = np.zeros((2*Lx*Ly, 2*Lx*Ly))

# Fill each energy level with up spin electron
for i in range(Ne):
  rhoE[i, i] = 1

print('rhoE is\n', rhoE)


# Construct H  = Ham_TB + Hsd
Ham_TB = hm.HTB(Lx, Ly, gam)
# print(np.kron(Ham_TB, np.eye(2)))


Ham_sd = hm.Hsd(Lx, Ly, pos_sp, Mag)
# print(Jsd*Ham_sd)

# Total Hamiltonian
Htot = np.kron(Ham_TB, np.eye(2)) + gfac*Jsd*Ham_sd


# Diagonalize the Hamiltonian
# Note in python, max eigvals are stored first
E, V = LA.eig(Htot)
# First eigenstate
# print(V[:, 0])

# Obtain rho in site basis : rho_site = V*rhoE*V^dag
# rho_site = np.dot(np.dot(V, rhoE), np.conjugate(V).T)
# rho_site = np.array([[1, 0], [0, 0]])
rho_site = rhoE 
# Note: if I convert to a different basis as Utkarsh mentioned, I don't
# get expected behavior of single z-electron interacting with x-spin 
# polarized classical spin

# print(rho_site)
# print('Trace is\n', np.trace(rho_site))

# Time evolution of rho(t) using rk4 
time = np.arange(0, tmax, dt)
# time = np.arange(0, 1)*dt

# S_opr
# Assumed that there are three classical spins  
sdim = Lx*Ly*2
sx_opr = np.zeros((sdim, sdim, sdim), dtype=complex)
sy_opr = np.zeros((sdim, sdim, sdim), dtype=complex)
sz_opr = np.zeros((sdim, sdim, sdim), dtype=complex)
 
for item in pos_sp:
  # print('item', item)
  proj = np.dot(hm.sitevec(Lx, Ly, item), hm.sitevec(Lx, Ly, item).T)
  sxop = np.kron(proj, SIG_X)
  syop = np.kron(proj, SIG_Y)
  szop = np.kron(proj, SIG_Z)

  sx_opr[item, :, :] = sxop
  sy_opr[item, :, :] = syop
  sz_opr[item, :, :] = szop

# Initial spin state at site 2
# Classical spin 

cspins1 = np.array(Mag[pos_sp[0]], dtype=float)
cspins2 = np.array(Mag[pos_sp[1]], dtype=float)
cspins3 = np.array(Mag[pos_sp[2]], dtype=float)
cspins4 = np.array(Mag[pos_sp[3]], dtype=float)
cspins = np.array(Mag, dtype = float)
# print(cspins[1])

# electron spin from rhoE initially
espins1 = np.array([0., 0., 0.], dtype=float)
espins2 = np.array([0., 0., 0.], dtype=float)
espins3 = np.array([0., 0., 0.], dtype=float)
espins4 = np.array([0., 0., 0.], dtype=float)

espins1[0] = np.trace(np.matmul(rho_site, sx_opr[0])).real
espins1[1] = np.trace(np.matmul(rho_site, sy_opr[0])).real
espins1[2] = np.trace(np.matmul(rho_site, sz_opr[0])).real

espins2[0] = np.trace(np.matmul(rho_site, sx_opr[1])).real
espins2[1] = np.trace(np.matmul(rho_site, sy_opr[1])).real
espins2[2] = np.trace(np.matmul(rho_site, sz_opr[1])).real

# TODO : Chane the line below when Jsd is added 

espins3[0] = np.trace(np.matmul(rho_site, sx_opr[2])).real
espins3[1] = np.trace(np.matmul(rho_site, sy_opr[2])).real
espins3[2] = np.trace(np.matmul(rho_site, sz_opr[2])).real

espins4[0] = np.trace(np.matmul(rho_site, sx_opr[3])).real
espins4[1] = np.trace(np.matmul(rho_site, sy_opr[3])).real
espins4[2] = np.trace(np.matmul(rho_site, sz_opr[3])).real

if method == 1:
 # rho at time t = 0
 rho_nw = rho_site
 # psi_e = np.array([1, 0]).reshape((2,1))  
 with open('../res/debug2.txt', "w") as f:
   for ind, t in enumerate(time):
     print('t = ', t)
     
     if run_LLG == 'yes':
       # Evolve classical spins using LLG solver

       B0x = B00 # *np.exp(-0.5*(t-t0)**2/Tw**2)
       cspins1 = llge.llg4(cspins1, espins1, cspins, t, dt, 
                                   Jsd, Jex, pos_sp[0], B0x, t0, Tw, B0z)
       cspins2 = llge.llg4(cspins2, espins2, cspins, t, dt, 
                                   Jsd, Jex, pos_sp[1], NoB0x, t0, Tw, -B0z)

       cspins3 = llge.llg4(cspins3, espins3, cspins, t, dt, 
                                   Jsd, Jex, pos_sp[2], NoB0x, t0, Tw, B0z)

       cspins4 = llge.llg4(cspins4, espins4, cspins, t, dt, 
                                   Jsd, Jex, pos_sp[3], NoB0x, t0, Tw, -B0z)
       Mag[pos_sp[0]] = list(cspins1)
       Mag[pos_sp[1]] = list(cspins2)
       Mag[pos_sp[2]] = list(cspins3)
       Mag[pos_sp[3]] = list(cspins4)

       cspins = np.array(Mag, dtype = float)

       # print('Mag', Mag[pos])
       Ham_sd = hm.Hsd(Lx, Ly, pos_sp, Mag)
       # print('Ham_sd is\n', gfac*Jsd*Ham_sd)
       # Total Hamiltonian
       Htot = np.kron(Ham_TB, np.eye(2)) + gfac*Jsd*Ham_sd
 

     # rk4 has bug : Fix the bug...seconded 
     rho_nw = rk4density.rk4(rho_nw, Htot, dt)
     
     # Update electron spin density at site 2
     espins1[0] = np.trace(np.matmul(rho_nw, sx_opr[0])).real
     espins1[1] = np.trace(np.matmul(rho_nw, sy_opr[0])).real
     espins1[2] = np.trace(np.matmul(rho_nw, sz_opr[0])).real

     espins2[0] = np.trace(np.matmul(rho_nw, sx_opr[1])).real
     espins2[1] = np.trace(np.matmul(rho_nw, sy_opr[1])).real
     espins2[2] = np.trace(np.matmul(rho_nw, sz_opr[1])).real

     espins3[0] = np.trace(np.matmul(rho_nw, sx_opr[2])).real
     espins3[1] = np.trace(np.matmul(rho_nw, sy_opr[2])).real
     espins3[2] = np.trace(np.matmul(rho_nw, sz_opr[2])).real

     espins4[0] = np.trace(np.matmul(rho_nw, sx_opr[3])).real
     espins4[1] = np.trace(np.matmul(rho_nw, sy_opr[3])).real
     espins4[2] = np.trace(np.matmul(rho_nw, sz_opr[3])).real

     # save time, cspins and espins  
     out_str = '%0.3f ' % t
     out_str +='%0.5e %0.5e %0.5e ' % (cspins1[0], cspins1[1], cspins1[2])
     out_str +='%0.5e %0.5e %0.5e ' % (cspins2[0], cspins2[1], cspins2[2])
     out_str +='%0.5e %0.5e %0.5e ' % (cspins3[0], cspins3[1], cspins3[2])
     out_str +='%0.5e %0.5e %0.5e ' % (cspins4[0], cspins4[1], cspins4[2])

     # Check that electron spin density is updating by plotting time evolution 
     out_str +='%0.5e %0.5e %0.5e ' % (espins1[0], espins1[1], espins1[2])
     out_str +='%0.5e %0.5e %0.5e ' % (espins2[0], espins2[1], espins2[2])
     out_str +='%0.5e %0.5e %0.5e ' % (espins3[0], espins3[1], espins3[2])
     out_str +='%0.5e %0.5e %0.5e ' % (espins4[0], espins4[1], espins4[2])
     out_str +='\n'
     f.write(out_str)

     # print('cspins2', cspins2)
