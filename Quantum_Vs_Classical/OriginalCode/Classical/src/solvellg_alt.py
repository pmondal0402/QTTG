import numpy as np
from numpy.linalg import norm as lanorm
from scipy.linalg import expm
from td2d.constants import SIG_X, SIG_Y, SIG_Z, HBAR
import Hamiltonians as hm
import rk4density 
from numpy import linalg as LA

gfac = 2 # Change to 2 

def generate_heff(Jsd, espins, time):
  # Landau g factor for spin

  heff_vals = []
  # Effective magnetic field due to sd coupling
  heff = np.array([0., 0., 0.], dtype=float)
  heff += -Jsd*gfac*espins
  # heff_vals.append(heff)
  return heff # np.array(heff_vals)

def dheun(cspins, espins, time, Jsd):
  hef = generate_heff(Jsd, espins, time)
  sh = np.cross(cspins, hef)
  ds = sh
  return sh

def norm_spins(spins):
  """ Normalizes spin vector"""
  norm_spins = lanorm(spins)
  nspins = spins/norm_spins
  return nspins

def llg4(cspins, espins, time, dt, Jsd):
  ds1 = dheun(cspins, espins, time, Jsd)*dt
  
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds1/2)

  ds2 = dheun(spins_t, espins, time_t, Jsd)*dt
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds2/2)

  ds3 = dheun(spins_t, espins, time_t, Jsd)*dt

  time_t = time + dt
  spins_t = norm_spins(cspins + ds3)

  ds4 = dheun(spins_t, espins, time_t, Jsd)*dt
  ds = (ds1 + 2*ds2 + 2*ds3 + ds4)/6.
  cspins = norm_spins(cspins + ds)  
  # cspins = norm_spins(cspins + 0.0*ds)
  return cspins

# Electronic time evolution
def generate_heffe(Jsd, espins, time):
  # Landau g factor for spin

  heff_vals = []
  # Effective magnetic field due to sd coupling
  heff = np.array([0., 0., 0.], dtype=float)
  heff += -Jsd*gfac*espins
  # heff_vals.append(heff)
  return heff # np.array(heff_vals)

def dheune(cspins, espins, time, Jsd):
  hef = generate_heffe(Jsd, espins, time)
  sh = np.cross(cspins, hef)
  ds = sh
  return sh

def llg4e(cspins, espins, time, dt, Jsd):
  ds1 = dheune(cspins, espins, time, Jsd)*dt
  
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds1/2)

  ds2 = dheune(spins_t, espins, time_t, Jsd)*dt
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds2/2)

  ds3 = dheune(spins_t, espins, time_t, Jsd)*dt

  time_t = time + dt
  spins_t = norm_spins(cspins + ds3)

  ds4 = dheune(spins_t, espins, time_t, Jsd)*dt
  ds = (ds1 + 2*ds2 + 2*ds3 + ds4)/6.
  cspins = norm_spins(cspins + ds)  
  # cspins = norm_spins(cspins + 0.0*ds)
  return cspins


# Overall plan TODO
# Params :
Lx = 1
Ly = 1

# Number of electronic spins
Ne = Lx*Ly 
# Hopping for tight binding 
gam = 1 # [eV]

# Params for sd Hamiltonian
# Number of classical spins
Nsp = 1 
# Jsd coupling 
Jsd = -0.9 

# Positions of classical spins
pos_sp = [0] # Choosing secind site 
# Vals for classical spins
# Mag = [[1, 0, 0], [1, 0, 0], [1, 0, 0]]
Mag = [[1, 0, 0]]

# Which method to use for time evolution
method = 1 # Note method 2 is working 

# rk4 params
dt = 0.1 # [fs]
tmax = 50  # [fs]

# Create rho in energy basis with half filled states
# i.e. |1up 2up 3up> for electrons
rhoE = np.zeros((2*Ne, 2*Ne))

# Fill each energy level with up spin electron
for i in range(0, 2*Ne, 2):
  rhoE[i, i] = 1

# Construct H  = Ham_TB + Hsd
Ham_TB = hm.HTB(Lx, Ly, gam)
# print(np.kron(Ham_TB, np.eye(2)))

Ham_sd = hm.Hsd(Lx, Ly, pos_sp, Mag)

# Todatl Hamiltonian
Htot = np.kron(Ham_TB, np.eye(2)) + Jsd*Ham_sd

# Diagonalize the Hamiltonian
# Note in python, max eigvals are stored first
E, V = LA.eig(Htot)
# First eigenstate
# print(V[:, 0])

# Obtain rho in site basis : rho_site = V*rhoE*V^dag
# rho_site = np.dot(np.dot(V, rhoE), np.conjugate(V).T)
rho_site = np.array([[1, 0], [0, 0]])
print(rho_site)
print('Trace is\n', np.trace(rho_site))

# TODO : Time evolution of rho(t) using rk4 
time = np.arange(0, tmax, dt)
# time = np.arange(0, 1)*dt

if method == 1:
   # rho at time t = 0
   rho_nw = rho_site
   # print('rho_old\n', rho_old)
   
   sx = [] ; sy = [] ; sz = []
   for ind, t in enumerate(time):
     # print(ind, t)
     rho_nw = rk4density.rk4(rho_nw, Htot, dt)
     # sx.append(np.trace(rho_nw*SIG_X).real)
     # sy.append(np.trace(rho_nw*SIG_Y).real)
     # sz.append(np.trace(rho_nw*SIG_Z).real)

    
     sx.append(np.trace(np.matmul(rho_nw, SIG_X)))
     sy.append(np.trace(np.matmul(rho_nw, SIG_Y)))
     sz.append(np.trace(np.matmul(rho_nw, SIG_Z)))
   # print('sx is\n', sx)
   print('sy is\n', sx) 
   # print('sz is\n', sz)

if method == 2:
   sx = [] ; sy = [] ; sz = []
   psi_e = np.array([1, 0]).reshape((2,1))
   for ind, t in enumerate(time):
      psi_t = np.matmul(expm(-1j*Ham_sd*t/HBAR), psi_e)
      psi_tc = np.conj(psi_t).T
      rho_nw = np.matmul(psi_t, psi_tc)
      
      # np.matmul(np.matmul(expm(1j*Htot*t/HBAR), rho_site),
      #   expm(-1j*Htot*t/HBAR))
      # sx.append(np.matmul(psi_tc,  np.matmul(SIG_X, psi_t))) 
      # sy.append(np.matmul(psi_tc,  np.matmul(SIG_Y, psi_t)) )
      # sz.append(np.matmul(psi_tc,  np.matmul(SIG_Z, psi_t)))

      sx.append(np.trace(np.matmul(rho_nw, SIG_X)))
      sy.append(np.trace(np.matmul(rho_nw, SIG_Y)))
      sz.append(np.trace(np.matmul(rho_nw, SIG_Z)))

   print('sy is\n', sy)
print('rho_new\n', rho_nw)

print('rho trace\n', np.matmul(rho_nw, SIG_Y))

print('rho trace\n', np.trace(np.matmul(rho_nw, SIG_Y)))
"""
# Example :
# Initial electronic state 
psi_e = np.array([1, 0]).reshape((2,1))
# psi_e  = psi_e/np.sqrt(2)

espins = np.array([0., 0., 1.], dtype=float)
# Evolve electron classically
classical_e = False

# Evolve electron quantum mechanically
quantum_e = True

Jsd = 0.1
dt = 0.1
time = 0.
# Initial spin state 
cspins = np.array([1., 0., 0.], dtype=float)

# Check if llg runs smoothly for fixed electronic spin 
tf = 50
times = np.arange(0, tf, dt)
with open('LLG_spin_elec.txt', "w") as f:
  for item in times:
    time = item
    # Evolve classical spins using LLG solver
    cspins_new = llg4(cspins, espins, time, dt, Jsd)
    cspins = cspins_new
    out_str = '%0.3f ' % time
    out_str +='%0.5e %0.5e %0.5e ' % (cspins[0], cspins[1], cspins[2])

    # Evolve electron as well
    # Solve Eq. 6 from Utkarsh's notes
    if classical_e :
      espins_new = llg4(espins, cspins, time, dt, Jsd)
      espins = espins_new
 
    # Alternate :
    elif quantum_e: 
    # Solve Schroedinger equation
    # Construct Hsd for electron
      Hsd = gfac*Jsd*( cspins[0]*SIG_X + cspins[1]*SIG_Y
                                     + cspins[2]*SIG_Z)
      # Electronic time evolution
      psi_t = np.matmul(expm(-1j*Hsd*time/HBAR), psi_e)
      psi_tc = np.conj(psi_t).T 
      # Measure and save electron spin density
      xx = np.matmul(psi_tc,  np.matmul(SIG_X, psi_t))
      yy = np.matmul(psi_tc,  np.matmul(SIG_Y, psi_t))
      zz = np.matmul(psi_tc,  np.matmul(SIG_Z, psi_t))

      espins[0] = xx[0][0]
      espins[1] = yy[0][0]
      espins[2] = zz[0][0]
    else:
     print('classical_e and quantum_e cannot be True (False) together')
 
    out_str +='%0.5e %0.5e %0.5e ' % (espins[0], espins[1], espins[2])
    out_str +='\n'
    f.write(out_str)
"""
