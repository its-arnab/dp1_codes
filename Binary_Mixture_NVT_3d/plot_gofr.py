# import os
# import glob
# import numpy as np
# import matplotlib.pyplot as plt

# # Base data path
# base_path = "../data_files"

# # Discover temperature folders
# temp_dirs = sorted(glob.glob(os.path.join(base_path, "T_*")))


# # Plot settings
# plt.figure(figsize=(8, 6))



# # Loop over temperature folders
# for temp_dir in temp_dirs:
#     # Extract temperature value from folder name
#     temp_str = os.path.basename(temp_dir).split("_")[1]
#     try:
#         temp_value = float(temp_str)
#     except ValueError:
#         continue  # skip malformed names

#     # Construct the expected MSD file path
#     gofr_file = os.path.join(temp_dir, f"gr_T_{temp_str}.dat")
#     if not os.path.isfile(gofr_file):
#         print(f"Missing: {gofr_file}")
#         continue

#     # Load and plot the data
#     try:
#         data = np.loadtxt(gofr_file)
#         distance = data[:, 0]
#         gofr = data[:, 1]
#         plt.plot(distance, gofr, linewidth=2, label=f"T = {temp_value:.2f}")
#         # plt.scatter(distance, gofr, s=10, alpha=1, label=f"T = {temp_value:.2f}")  # Scatter plot for every 10th point


#     except Exception as e:
#         print(f"Error reading {gofr_file}: {e}")
#         continue

# # Finalize plot

# plt.xlabel("Distance", fontsize=16)
# plt.ylabel(r"$g(r)$", fontsize=16)
# plt.legend(fontsize=12)
# plt.grid(True, which='major', linestyle='--', linewidth=0.4)
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.tight_layout()
# plt.show()


import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Base data path
base_path = "../data_files"
temp_dirs = sorted(glob.glob(os.path.join(base_path, "T_*")))

# Main plot setup
fig, ax = plt.subplots(figsize=(8, 6))

# Create inset axes
ax_inset1 = inset_axes(ax, width="30%", height="30%", loc='upper center')
                    #    bbox_to_anchor=(1.75, 2.45, 1, 1), bbox_transform=ax.transAxes)
ax_inset2 = inset_axes(ax, width="30%", height="30%", loc='center right')
                    #    bbox_to_anchor=(1.75, 1.02, 1, 1), bbox_transform=ax.transAxes)

# Loop through temperature directories
for temp_dir in temp_dirs:
    temp_str = os.path.basename(temp_dir).split("_")[1]
    try:
        temp_value = float(temp_str)
    except ValueError:
        continue

    gofr_file = os.path.join(temp_dir, f"gr_T_{temp_str}.dat")
    if not os.path.isfile(gofr_file):
        print(f"Missing: {gofr_file}")
        continue

    try:
        data = np.loadtxt(gofr_file)
        distance = data[:, 0]
        gofr = data[:, 1]
        ax.plot(distance, gofr, linewidth=2, label=f"T = {temp_value:.2f}")
        ax_inset1.plot(distance, gofr, linewidth=1)
        ax_inset2.plot(distance, gofr, linewidth=1)
    except Exception as e:
        print(f"Error reading {gofr_file}: {e}")
        continue

# Main plot formatting
ax.set_xlabel("Distance(r)", fontsize=16)
ax.set_ylabel(r"$g(r)$", fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, which='major', linestyle='--', linewidth=0.4)
ax.tick_params(axis='both', which='major', labelsize=14)

# Inset formatting
ax_inset1.set_xlim(1.0, 1.15)
ax_inset1.set_ylim(2.4, 3.30)
ax_inset1.grid(True, which='major', linestyle='--', linewidth=0.4)
ax_inset1.tick_params(axis='both', which='both', labelsize=8, direction='in', length=3)


ax_inset2.set_xlim(1.2, 1.7)
ax_inset2.set_ylim(0.2, 0.7)
ax_inset2.axvline(x=1.432, ymax=0.7, c='k', ls='--', lw=1.0, label='$r_{caged}=1.432$')
ax_inset2.grid(True, which='major', linestyle='--', linewidth=0.4)
ax_inset2.tick_params(axis='both', which='both', labelsize=8, direction='in', length=3)
ax_inset2.legend(loc='upper center', fontsize=10)

plt.show()
