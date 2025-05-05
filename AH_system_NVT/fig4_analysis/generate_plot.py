import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterSciNotation

# Get all output files
files = glob.glob("avg_data/N_*_rho_*_K_*_T_*_dt_*")

# Data dictionary to store energy and pressure separately
data = {}

for file in files:
    # Extract parameters from filename
    parts = file.split("_")
    N = int(parts[2])   # Extract N value
    rho = float(parts[4]) 
    K = float(parts[6]) # Extract K value
    T = float(parts[8]) # Extract T value

    # Read both avg energy and avg pressure from the file
    with open(file, "r") as f:
        lines = f.readlines()
        avg_energy = float(lines[0].strip())  # First line: Energy
        avg_pressure = float(lines[1].strip())  # Second line: Pressure

    # Store data separately for Energy and Pressure
    if (N, K) not in data:
        data[(N, K)] = {"T": [], "E": [], "P": []}  # Create empty lists

    data[(N, K)]["T"].append(T)
    data[(N, K)]["E"].append(avg_energy)
    data[(N, K)]["P"].append(avg_pressure)

# Plot Energy vs. Temperature for different N, K
plt.figure(1, figsize=(8, 8))
for (N, K), values in data.items():
    sorted_indices = np.argsort(values["T"])  # Sort T values for smooth plotting
    T_sorted = np.array(values["T"])[sorted_indices]
    E_sorted = np.array(values["E"])[sorted_indices]

    plt.semilogx(T_sorted, E_sorted, marker='o', label=f"N={N}, K={K}")

# Generate exactly 5 log-spaced x-ticks
x_ticks = np.array([0.01, 0.1, 1.0])
plt.xticks(x_ticks)
plt.plot([0.4, 0.4], [-2.8, 2.0], lw = 0.8, ls='--', c='k')
plt.gca().set_box_aspect(1)
plt.xlabel("Temperature (T)")
plt.ylabel(r"$ \frac{\langle E \rangle} {N}$")
plt.legend()
plt.title("Energy vs. Temperature")

# Plot Pressure vs. Temperature for different N, K
plt.figure(2, figsize=(8, 8))
for (N, K), values in data.items():
    sorted_indices = np.argsort(values["T"])
    T_sorted = np.array(values["T"])[sorted_indices]
    P_sorted = np.array(values["P"])[sorted_indices]
    valid_indices = P_sorted > 0
    plt.loglog(T_sorted[valid_indices], P_sorted[valid_indices], marker='s', label=f"N={N}, K={K}")

# Set x-ticks on a log scale
x_ticks = np.array([0.01, 0.1, 1.0])
plt.xticks(x_ticks)
y_ticks = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1])
plt.yticks(y_ticks)
plt.plot([0.4, 0.4], [1e-4, 0.1], lw = 0.8, ls ='--', c='k')
plt.gca().set_box_aspect(1)
plt.xlabel("Temperature (T)")
plt.ylabel("Average Pressure (P)")
plt.legend()
plt.title("Pressure vs. Temperature")
plt.show()
