import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

processors = [1, 2, 4]
times = np.genfromtxt("times.dat").T
meaning = ["Real", "User", "System"]

fig, axes = plt.subplots(3, sharex=True)
for i in range(3):
    axes[i].plot(processors, times[i], "-o", label = meaning[i])
    axes[i].set_ylabel("Time (s)")
    axes[i].legend()
    axes[i].get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))
axes[-1].set_xlabel("Processors")
plt.savefig("time.pdf")
