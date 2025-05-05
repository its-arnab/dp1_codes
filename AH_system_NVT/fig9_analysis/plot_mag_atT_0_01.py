import numpy as np
import matplotlib.pyplot as plt

N = np.array([200, 1000, 2000, 4000, 6000, 8000, 10000 ], dtype='float')
#N = np.log10(N)
mag = np.array([0.996513, 0.993705, 0.917701, 0.567846, 0.44639, 0.3532, 0.32497])

plt.figure(figsize=(4,4))
plt.plot(N, mag, 'o-', c='blue', lw = 1.5, markersize=4)
plt.xlabel(r'$N$')
plt.ylabel(r'$ \langle M \rangle $')
plt.grid(which='major')
plt.ylim(0.2, 1.05)
plt.savefig("mag_plot_T_0_01.pdf", bbox_inches='tight', dpi=150)
plt.show()
