import numpy as np
from numpy.linalg import norm as lanorm
from scipy.linalg import expm
from td2d.constants import SIG_X, SIG_Y, SIG_Z, HBAR

gfac = 2 # Change to 2 

def generate_heff(Jsd, espins, time):
  # Landau g factor for spin

  heff_vals = []
  # Effective magnetic field due to sd coupling
  heff = np.array([0., 0., 0.], dtype=float)
  heff += -Jsd*gfac*espins
  heff_vals.append(heff)
  return np.array(heff_vals)

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
  # cspins = norm_spins(cspins + ds)  
  cspins = norm_spins(cspins + 0.0*ds)
  return cspins


# Example :
# Initial electronic state 
psi_e = np.array([1, 0]).reshape((2,1))
# psi_e  = psi_e/np.sqrt(2)

espins = np.array([0., 0., 1.], dtype=float)
Jsd = 0.05/2
dt = 0.01
time = 0.
# Initial spin state 
cspins = np.array([1., 0., 0.], dtype=float)

# Check if llg runs smoothly for fixed electronic spin 
tf = 500
times = np.arange(0, tf, dt)
with open('LLG_spin_elec.txt', "w") as f:
  for item in times:
    time = item
    # Evolve classical spins using LLG solver
    cspins_new = llg4(cspins, espins, time, dt, Jsd)
    cspins = cspins_new
    out_str = '%0.3f ' % time
    out_str +='%0.5e %0.5e %0.5e ' % (cspins[0, 0], cspins[0, 1], cspins[0, 2])

    # Evolve electron as well
    # Construct Hsd for electron
    Hsd = gfac*Jsd*( cspins[0, 0]*SIG_X + cspins[0, 1]*SIG_Y
                                   + cspins[0, 2]*SIG_Z)
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
 
    out_str +='%0.5e %0.5e %0.5e ' % (espins[0], espins[1], espins[2])
    out_str +='\n'
    f.write(out_str)
