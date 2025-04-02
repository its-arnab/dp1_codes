import os
import numpy as np
import matplotlib.pyplot as plt
import glob

# Folder containing the files
folder_path = "Vcm_data/"  # Change this to your folder path

# Get list of files matching the pattern T_%.2f
file_list = sorted(glob.glob(os.path.join(folder_path, "T_*")))  # Change extension if needed

temperatures = []
all_vx_data = {}

for file in file_list:
    # Extract temperature from filename
    try:
        temp = float(os.path.basename(file).split('_')[1])
        temperatures.append(temp)
    except ValueError:
        print(f"Skipping file {file}, unable to extract temperature.")
        continue
    
    # Read the first column (vx) from the file
    data = np.loadtxt(file)
    vx = data[:, 0]  # First column
    all_vx_data[temp] = vx

# Sort temperatures for ordered plotting
temperatures.sort()

# Plot histograms for each temperature
plt.figure(figsize=(8, 6))

for temp in temperatures:
    vx = all_vx_data[temp]
    
    # Compute histogram
    counts, bin_edges = np.histogram(vx, bins=50, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Plot in semi-log scale
    plt.semilogy(bin_centers, counts, label=f"T = {temp:.2f}")

plt.xlabel(r"$v_x$", fontsize=14)
plt.ylabel(r"$P(v_x)$", fontsize=14)
plt.title("Histogram of vx for Different Temperatures")
plt.legend()
plt.grid()
plt.show()