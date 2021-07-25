from td2d.constants import HBAR
from numpy import matmul

def rk4(rho, H, dt):
  # Runge-Kutta method to evolve :
  # rho'(t) = (1i/hbar)*[rho(t),H];
  fac = 1j/HBAR
  y0 = rho
  k1 =  fac*commutation(y0, H)
  y1 = rho + 0.5*dt*k1
  
  k2 = fac*commutation(y1, H) 
  y2 = rho + 0.5*dt*k2

  k3 = fac*commutation(y2, H) 
  y3 = rho + dt*k3

  k4 = fac*commutation(y3, H) 
  # print((dt/6.0)*(k1 + 2*k2 + 2*k3 + k4))
  res = rho + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
  return res
 

def commutation(A, B):
  # Returns matrix commutation
  res = matmul(A, B) - matmul(B, A)
  return res
