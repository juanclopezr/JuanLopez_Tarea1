import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#data = np.genfromtxt('atpos.dat').T
data = np.loadtxt(Energies.dat)

plt.plot(data[:,0],data[:,1])
plt.plot(data[:,0],data[:,2])
plt.plot(data[:,0],data[:,1])
plt.xlabel("time")
plt.ylabel("Energy")

#plt.imshow(data, cmap='jet')
#plt.xlabel("time")
#plt.ylabel("particle")
#plt.set_cmap('nipy_spectral')
plt.savefig('graph.pdf')
