import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

# Constants
V  = 1.0e-15	# L
NA = 6.022e23	# molecules/mole
tstart =  0.0	# s
tend   = 30.0	# s

# Rates
kf=1.07e5/(NA*V) # /M/s
kr=0.351	 # /s

# Initial Species Counts
A  = 1000
B  = 1000
C  = 0
S0 = [A, B, C]

# Definition of ODEs
def ds_dt(s, t):
	Ai = s[0]
	Bi = s[1]
	Ci = s[2]
	# Rate equations
	dA_dt = -kf*Ai*Bi + kr*Ci
	dB_dt = -kf*Ai*Bi + kr*Ci
	dC_dt =  kf*Ai*Bi - kr*Ci
	return [dA_dt, dB_dt, dC_dt]

# Solve
t    = np.linspace(tstart, tend, 1000000)
soln = spi.odeint(ds_dt, S0, t)

# Plot
plt.figure()
plt.plot(t, soln[:,0], label="A/B")
plt.plot(t, soln[:,2], label="C")
plt.xlabel('Time (s)')
plt.ylabel('Molecule Count')
plt.legend()
plt.savefig('BimolecularODE.png')


