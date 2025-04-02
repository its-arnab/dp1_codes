import numpy as np
import matplotlib.pyplot as plt
import re
import os

def plot_particles_from_file(filename):
    # Extract parameters from filename
    match = re.search(r'N_(\d+)_T_(\d+\.\d+)', os.path.basename(filename))
    if match:
        N = int(match.group(1))
        T = float(match.group(2))
    else:
        N = 0
        T = 0
        print("Could not extract parameters from filename")
    
    # Read the data file
    try:
        data = np.loadtxt(filename)
        
        # Extract columns
        x = data[:, 0]
        y = data[:, 1]
        theta = data[:, 2]
        
        # Calculate L based on data (assuming square domain)
        L = max(np.max(x), np.max(y))
        
        # Create the plot
        color_values = np.cos(theta)
        fig, ax = plt.subplots(figsize=(8, 8))
        sc = ax.scatter(x, y, c=color_values, cmap='RdBu', vmin=-1, vmax=1, s=10)
        ax.set_xlim(0, L+1)
        ax.set_ylim(0, L+1)
        cbar = plt.colorbar(sc)
        cbar.set_label(r'$\cos(\theta_i)$')
        plt.gca().set_box_aspect(1)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title(f"AH Model Evolution (N={N}, T={T})")
        
        # Show the plot
        plt.show()
        
    except Exception as e:
        print(f"Error reading or plotting data: {e}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        # If a filename is provided as command line argument
        filename = sys.argv[1]
    else:
        # Default filename
        filename = 'final_config/N_00200_T_0.40357'
    
    plot_particles_from_file(filename)
