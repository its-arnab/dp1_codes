import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Base data path
base_path = "../data_files"

# Discover temperature folders
temp_dirs = sorted(glob.glob(os.path.join(base_path, "T_*")))

# Plot settings
plt.figure(figsize=(8, 6), dpi=100)

# Loop over temperature folders
for temp_dir in temp_dirs:
    # Extract temperature value from folder name
    temp_str = os.path.basename(temp_dir).split("_")[1]
    try:
        temp_value = float(temp_str)
    except ValueError:
        continue  # skip malformed names
    
    # Construct the expected MSD file path
    msd_file = os.path.join(temp_dir, f"msd_T_{temp_str}.dat")
    if not os.path.isfile(msd_file):
        print(f"Missing: {msd_file}")
        continue
    
    # Load and plot the data
    try:
        data = np.loadtxt(msd_file)
        times = data[:, 0]
        msd = data[:, 1]
        plt.plot(times, msd, linewidth=2, label=f"T = {temp_value:.3f}")
    except Exception as e:
        print(f"Error reading {msd_file}: {e}")
        continue

# Add reference lines
# For ballistic region (t²) - selecting short time range
t_ballistic = np.logspace(-2.5, -0.5, 100)  # Adjust range as needed
# Scale the t² line appropriately to overlay with your data
ballistic_scale = 0.8  # Adjust this scaling factor to fit your data
msd_ballistic = ballistic_scale * t_ballistic**2
plt.plot(t_ballistic, msd_ballistic, 'k--', linewidth=2, label=r'$\sim t^2$ (Ballistic)')

# For diffusive region (t¹) - selecting longer time range
t_diffusive = np.logspace(1, 3, 100)  # Adjust range as needed
# Scale the t¹ line appropriately to overlay with your data
diffusive_scale = 0.1  # Adjust this scaling factor to fit your data
msd_diffusive = diffusive_scale * t_diffusive
plt.plot(t_diffusive, msd_diffusive, 'k-.', linewidth=2, label=r'$\sim t$ (Diffusive)')

# Finalize plot
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Time", fontsize=16)
plt.ylabel("MSD", fontsize=16)
plt.xlim(0.005, 10000)
plt.ylim(1e-5, 1e2)
plt.legend(fontsize=12)
plt.grid(True, which='major', linestyle='--', linewidth=0.5)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.show()
