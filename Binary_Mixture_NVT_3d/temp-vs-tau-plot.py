'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Known constants
tau0 = 0.3101
K_vft = 0.2243

# Define the VFT function with only T_vft as variable
def vft_fit(T, T_vft):
    return tau0 * np.exp(1.0 / (K_vft * (T / T_vft - 1.0)))

# Data
tau_alpha = np.array([1.65, 2.16, 2.96, 4.63, 8.69, 24.91, 229, 3223])
tau_alpha2 = np.array([0.99, 1.24, 1.67, 2.61, 4.91, 13.51, 125.33, 1908.5])
temperature = np.array([1.1, 1.0, 0.90, 0.80, 0.70, 0.60, 0.50, 0.44])

# Fit for tau_alpha
popt1, _ = curve_fit(vft_fit, temperature, tau_alpha, p0=[0.3])
T_vft_fit1 = popt1[0]

# Fit for tau_alpha2
popt2, _ = curve_fit(vft_fit, temperature, tau_alpha2, p0=[0.3])
T_vft_fit2 = popt2[0]

# Generate fitted curves
T_fine = np.linspace(0.42, 1.12, 300)
fit_curve1 = vft_fit(T_fine, T_vft_fit1)
fit_curve2 = vft_fit(T_fine, T_vft_fit2)

# Plot
plt.figure(figsize=(6,5))
plt.semilogy(temperature, tau_alpha, 'o', color='red', label=r'$\tau_{\alpha}$ data')
plt.semilogy(temperature, tau_alpha2, '^', color='blue', label=r'From $F_s(k,t)$')

plt.semilogy(T_fine, fit_curve1, '--', color='red', label=fr'VFT fit 1: $T_{{\rm vft}}$ = {T_vft_fit1:.3f}')
plt.semilogy(T_fine, fit_curve2, '--', color='blue', label=fr'VFT fit 2: $T_{{\rm vft}}$ = {T_vft_fit2:.3f}')

plt.ylabel(r'$\tau_{\alpha}$ ', fontsize=14)
plt.xlabel('Temperature', fontsize=14)
plt.tick_params(labelsize=14)
plt.grid(True, alpha=0.3, which='major')
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

# Print fitted parameters
print(f"Fitted T_vft for tau_alpha: {T_vft_fit1:.4f}")
print(f"Fitted T_vft for tau_alpha2: {T_vft_fit2:.4f}")
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the VFT function
def vft_full(T, tau0, K_vft, T_vft):
    return tau0 * np.exp(1.0 / (K_vft * (T / T_vft - 1.0)))

# Data
tau_alpha = np.array([1.65, 2.16, 2.96, 4.63, 8.69, 24.91, 229, 3223])
tau_alpha2 = np.array([0.99, 1.24, 1.67, 2.61, 4.91, 13.51, 125.33, 1908.5])
temperature = np.array([1.1, 1.0, 0.90, 0.80, 0.70, 0.60, 0.50, 0.44])

# Initial guesses and bounds: [tau0, K_vft, T_vft]
p0 = [0.3, 0.2, 0.3]
bounds = ([1e-4, 0.01, 0.1], [10.0, 1.0, 0.6])  # reasonable physical bounds

# Fit for tau_alpha
popt1, _ = curve_fit(vft_full, temperature, tau_alpha, p0=p0, bounds=bounds)
tau0_1, K_1, T_vft_1 = popt1

# Fit for tau_alpha2
popt2, _ = curve_fit(vft_full, temperature, tau_alpha2, p0=p0, bounds=bounds)
tau0_2, K_2, T_vft_2 = popt2

# Generate fitted curves
T_fine = np.linspace(0.42, 1.12, 300)
fit_curve1 = vft_full(T_fine, *popt1)
fit_curve2 = vft_full(T_fine, *popt2)

# Plot
plt.figure(figsize=(6,5))
plt.semilogy(temperature, tau_alpha, 'o', color='red', label=r'$\tau_{\alpha}$ from $Q(t)$')
plt.semilogy(temperature, tau_alpha2, '^', color='blue', label=r'$\tau_{\alpha}$ from $F_s(k,t)$')

plt.semilogy(T_fine, fit_curve1, '--', color='red', label=fr'Fit 1: $T_{{\rm vft}}$={T_vft_1:.3f}, $K$={K_1:.3f}')
plt.semilogy(T_fine, fit_curve2, '--', color='blue', label=fr'Fit 2: $T_{{\rm vft}}$={T_vft_2:.3f}, $K$={K_2:.3f}')

plt.ylabel(r'$\tau_{\alpha}$ ', fontsize=14)
plt.xlabel('Temperature', fontsize=14)
plt.tick_params(labelsize=14)
plt.grid(True, alpha=0.3, which='major')
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()

# Print fitted parameters
print("Fitted parameters:")
print(f"From Q(t): tau0 = {tau0_1:.4e}, K = {K_1:.4f}, T_vft = {T_vft_1:.4f}")
print(f"From Fs(k,t): tau0 = {tau0_2:.4e}, K = {K_2:.4f}, T_vft = {T_vft_2:.4f}")


