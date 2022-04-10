import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc


Nx = 101;                         # Number of grid points
xmax = 2.;                        # Domain limit to the right
xmin = -2.;                       # Domain limit to the left
dx = (xmax-xmin)/(Nx-1)           # Mesh size
x = np.linspace(xmin,xmax,Nx)     # Discretized mesh

dt = 0.04                         # Time step
t_end = 5.                        # Final time
Nt = int(t_end/dt)                # Number of iterations
t = np.linspace(0.,t_end,Nt+1)    # Time vector


c = 0.8                           # Advection speed
CFL = c*dt/dx                     # CFL number
U = np.zeros((Nt+1,Nx))           # u^n_i
U[0,:] = np.exp(-0.5*(x/0.4)**2)  # Initial solution
Uex = U[0,:]                      # Exact solution

#===============================================================
# Compute exact solution
#===============================================================
for n in range (0,Nt):
  d = c*(n+1)*dt
  Uex = np.exp(-0.5*(np.mod(x-d+xmax,4)-xmax)**2/0.4**2)
  errL1 = U - Uex
  errL2 = np.linalg.norm(errL1)

#===============================================================
# Plot solution
#===============================================================

  if (n==0): fig, ax = plt.subplots(figsize=(5.5,4))
  plt.clf()
  plt.scatter(x,Uex, marker='o', facecolors='white', color='k')
  plt.legend('Exact solution')
  plt.axis([xmin, xmax, 0, 1.4])
  plt.title('t='+str(round(dt*(n+1),3)),fontsize=16)
  plt.xlabel('x',fontsize=18)
  plt.ylabel('u',fontsize=18)
  plt.subplots_adjust(left=0.2)
  plt.subplots_adjust(bottom=0.18)
  plt.draw()
  plt.pause(0.001)
plt.show()