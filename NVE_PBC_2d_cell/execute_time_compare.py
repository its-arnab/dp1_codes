import numpy as np
import matplotlib.pyplot as plt

# Given Data
N_values = np.array([100, 196, 400, 576, 900, 1225, 1600, 2025, 2500, 3600, 4900, 6400])
N2_algo      = np.array([0.055, 0.153, 0.518, 0.969, 2.092, 3.505, 5.625, 8.587, 12.888, 29.292, 62.379, 101.335])
verlet_list  = np.array([0.054, 0.085, 0.186, 0.343, 0.492, 0.710, 1.001, 1.313, 1.779, 3.159, 4.802, 7.507])
cell_list_naive = np.array([0.044, 0.091, 0.222, 0.288, 0.447, 0.700, 0.975, 1.213, 1.398, 1.922, 2.882, 3.532])
cell_list_pro   = np.array([0.046, 0.092, 0.199, 0.264, 0.450, 0.579, 0.779, 0.956, 1.218, 1.807, 2.511, 3.301])

# Create a mask to filter data for N in [100, 1000]
mask = (N_values >= 100) & (N_values <= 1000)

# Create two subplots (side by side)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=False)

# -------------------------------
# Full-range plot (Left subplot)
# -------------------------------
ax1.plot(N_values, N2_algo, 'o-', label="N² Algorithm", linewidth=2)
ax1.plot(N_values, verlet_list, 's-', label="Verlet List", linewidth=2)
# ax1.plot(N_values, cell_list_naive, 'd-', label="Cell List Naive", linewidth=2)
ax1.plot(N_values, cell_list_pro, '^-', label="Cell List", linewidth=2)

ax1.set_xlim(0, 6500)
ax1.set_ylim(0, 10)  # Adjust y-limit if needed for full-range visualization
ax1.set_xlabel("Nuber of Particles (N)", fontsize=12)
ax1.set_ylabel("Execution Time (s)", fontsize=12)
ax1.set_title("Full-range Execution Time", fontsize=14)
ax1.grid(True, which="both", linestyle="--", linewidth=0.6)
ax1.legend(fontsize=10)

# -------------------------------
# Zoomed plot (Right subplot)
# -------------------------------
ax2.plot(N_values, N2_algo, 'o-', label="N² Algorithm", linewidth=2)
ax2.plot(N_values, verlet_list, 's-', label="Verlet List", linewidth=2)
# ax2.plot(N_values, cell_list_naive, 'd-', label="Cell List Naive", linewidth=2)
ax2.plot(N_values, cell_list_pro, '^-', label="Cell List", linewidth=2)

# Focus on the desired region: N=100 to N=1000
ax2.set_xlim(0, 600)
# Set y-limit based on the zoomed data range; adjust as needed
ax2.set_ylim(0, 0.5)  

ax2.set_xlabel("Number of Particles (N)", fontsize=12)
ax2.set_title("Zoomed Region (N = 100 to 1000)", fontsize=14)
ax2.grid(True, which="both", linestyle="--", linewidth=0.6)
ax2.legend(fontsize=10)
# Improve layout spacing
plt.tight_layout()
plt.show()
