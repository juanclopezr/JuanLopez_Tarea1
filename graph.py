import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('atpos.dat')

print(data.shape)
plt.imshow(data, cmap='jet')
plt.show()
