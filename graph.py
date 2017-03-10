import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t, e1, e2, e3 = np.genfromtxt("Energies.dat").T
plt.plot(t, e1,label="$E_1$")
plt.plot(t, e2, label="$E_2$")
plt.plot(t, e3, label="$E_3$")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig('graph.pdf')
