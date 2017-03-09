import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = np.genfromtxt("Energies.dat")

plt.plot(data[:,0],data[:,1], "-", ms=1, label="$E_1$")
plt.plot(data[:,0],data[:,2], "-", ms=1, label="$E_2$")
plt.plot(data[:,0],data[:,3], "-", ms=1, label="$E_3$")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig('graph.pdf')
