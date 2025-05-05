import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Base data path
base_path = "../data_files"

# Discover temperature folders
temp_dirs = sorted(glob.glob(os.path.join(base_path, "T_*")))


# Plot settings
plt.figure(figsize=(8, 6))


# find tau_alpha
tau_alpha =[]
def find_tau_alpha(data):
    # Find the index where Q(t) drops below 1/e
    index = np.where(data[:, 1] < 1/np.e)[0]
    if index.size > 0:
        return data[index[0], 0]  # Return the time value at that index
    else:
        return None

# Loop over temperature folders
for temp_dir in temp_dirs:
    # Extract temperature value from folder name
    temp_str = os.path.basename(temp_dir).split("_")[1]
    try:
        temp_value = float(temp_str)
    except ValueError:
        continue  # skip malformed names

    # Construct the expected MSD file path
    fskt_file = os.path.join(temp_dir, f"fskt_k00725_T_{temp_str}.dat")
    if not os.path.isfile(fskt_file):
        print(f"Missing: {fskt_file}")
        continue

    # Load and plot the data
    try:
        data = np.loadtxt(fskt_file)
        times = data[:, 0]
        qoft = data[:, 1]
        plt.plot(times, qoft, linewidth=2, label=f"T = {temp_value:.2f}")
        # plt.scatter(times, qoft, s=10, alpha=1, label=f"T = {temp_value:.2f}")  # Scatter plot for every 10th point

        # store tau_alpha
        tau = find_tau_alpha(data)
        if tau is not None:
            tau_alpha.append(tau)
            print(f"tau_alpha for T = {temp_value:.2f} is {tau:.4f}")
        else:
            print(f"tau_alpha not found for T = {temp_value:.2f}")

    except Exception as e:
        print(f"Error reading {fskt_file}: {e}")
        continue

# Finalize plot
plt.xscale("log")
plt.xlabel("Time", fontsize=16)
plt.ylabel(r"$F_s(k,t)$", fontsize=16)
# plt.xlim(0.005, 10000)
plt.axhline(y=1/np.e, color='k', linestyle='--', label='1/e')
plt.legend(fontsize=12)
plt.grid(True, which='major', linestyle='--', linewidth=0.4)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.show()





'''
At higher temperatures, the structural relaxation time τ_α is smaller because particles have more thermal energy (~k_BT),
enabling them to move and rearrange more freely. This increased kinetic energy helps them overcome local energy barriers, 
break out of cages formed by neighboring particles, and relax the structure quickly. 

As a result, the system loses its memory of initial configurations faster, leading to a rapid decay of Q(t). 
In contrast, at low temperatures, motion is restricted, cage-breaking is rare, and relaxation becomes sluggish, causing τ_α to grow significantly.
'''
