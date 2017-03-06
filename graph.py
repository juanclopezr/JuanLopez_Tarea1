import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = np.genfromtxt('atpos.dat').T

plt.imshow(data, cmap='jet')
plt.xlabel("time")
plt.ylabel("particle")
#plt.set_cmap('nipy_spectral')
plt.savefig('graph.pdf')
