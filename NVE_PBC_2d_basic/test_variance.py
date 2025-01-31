import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt(r"dt_0.001_data.txt")
data2 = np.loadtxt(r"dt_0.002_data.txt")

x1 = data1[:,0]*0.001
x2 = data2[:,0]*0.002

y1 = data1[:,4]
y2 = data2[:,4]

del data1
del data2


plt.plot(x1, y1, lw = 1.0, label='dt=0.001')
plt.plot(x2, y2, lw = 1.0, label='dt=0.002')

plt.axhline(y = np.mean(y1), ls = '--', lw = 0.8)
plt.axhline(y = np.mean(y2), ls = '--', lw = 0.8)

temp = np.var(y2)/np.var(y1)
print('variance ratio for two case:', np.sqrt(temp))

plt.legend();
plt.show()

