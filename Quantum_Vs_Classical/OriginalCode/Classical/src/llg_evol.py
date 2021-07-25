import numpy as np
from numpy.linalg import norm as lanorm

gfac = 2 # Change to 2 

def generate_heff(Jsd, Jex, espins, cspinstot, time, ind, B0x, t0, Tw, B0z):
  # Landau g factor for spin

  # TODO : Add pulse Hamiltonian along x-dir on Spin 1
  # Check if this vector is along x-dir
  xvec = 0*cspinstot[0]
  yvec = 0*cspinstot[0]
  zvec = 0*cspinstot[0]
  xvec[0] = 1
  yvec[1] = 1
  zvec[2] = 1 

  heff_vals = []

  # Effective magnetic field due to sd coupling
  heff = np.array([0., 0., 0.], dtype=float)
  heff += -Jsd*gfac*espins
  # heff_vals.append(heff)
  # print('ind is ', ind)
  # Add exchange coupling 
  # Find nearest neighbour

  if ind == 0:
    # print('Do nothing ', cspinstot)
    # nearest neighbours
    heff+= -Jex* cspinstot[ind+1]

    heff += B0z*zvec
    # Add external Bx only at site 0
    heff +=B0x*xvec

  if ind == 1:
    # nearest neighbours
    heff+= -Jex*(cspinstot[ind-1] + cspinstot[ind+1])

    heff += B0z*zvec

    # Add external Bx only at site 0
    heff +=B0x*xvec


  if ind == 2:
    # nearest neighbours
    heff+= -Jex*(cspinstot[ind-1] + cspinstot[ind+1])

    heff += B0z*zvec

    # Add external Bx only at site 0
    heff +=B0x*xvec

  if ind == 3:
    # print('Do nothing ', cspinstot)
    # nearest neighbours
    heff+= -Jex*cspinstot[ind-1]
    heff += B0z*zvec

  return heff # np.array(heff_vals)

def dheun(cspins, espins, cspinstot, time, Jsd, Jex, ind, B0x, t0, Tw, B0z):
  hef = generate_heff(Jsd, Jex, espins, cspinstot, time, ind, B0x, t0, Tw, B0z)
  sh = np.cross(cspins, hef)
  ds = sh
  return sh

def norm_spins(spins):
  """ Normalizes spin vector"""
  norm_spins = lanorm(spins)
  nspins = spins/norm_spins
  return nspins

def llg4(cspins, espins, cspinstot, time, dt, Jsd, Jex, ind, B0x, t0, Tw, B0z):
  ds1 = dheun(cspins, espins, cspinstot, time, Jsd, Jex, ind, B0x, t0, Tw, B0z)*dt
  
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds1/2)

  ds2 = dheun(spins_t, espins, cspinstot, time_t, Jsd, Jex, ind, B0x, t0, Tw,
                                                                      B0z)*dt
  time_t = time + dt/2
  spins_t = norm_spins(cspins + ds2/2)

  ds3 = dheun(spins_t, espins, cspinstot, time_t, Jsd, Jex, ind, B0x, t0, Tw,
                                                                      B0z)*dt

  time_t = time + dt
  spins_t = norm_spins(cspins + ds3)

  ds4 = dheun(spins_t, espins, cspinstot, time_t, Jsd, Jex, ind, B0x, t0, Tw,
                                                                      B0z)*dt
  ds = (ds1 + 2*ds2 + 2*ds3 + ds4)/6.
  cspins = norm_spins(cspins + ds)  
  # cspins = norm_spins(cspins + 0.0*ds)
  return cspins

