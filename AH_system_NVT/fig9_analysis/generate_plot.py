import os
import numpy as np
import matplotlib.pyplot as plt

def read_magnetisation_data(folder="magnetisation_data"):
    data_dict = {}
    
    # Check if the folder exists
    if not os.path.exists(folder):
        print(f"Error: Folder '{folder}' not found.")
        return data_dict
    
    # Look for folders with pattern N_XXXXX
    for item in os.listdir(folder):
        item_path = os.path.join(folder, item)
        
        # Check if it's a directory and matches the pattern N_xxxxx
        if os.path.isdir(item_path) and item.startswith("N_"):
            try:
                # Extract N value from the folder name
                N_value = int(item[2:])
                
                # Look for mag_data file in this directory
                mag_data_path = os.path.join(item_path, "mag_data")
                
                if os.path.exists(mag_data_path):
                    # Read the mag_data file
                    data = np.loadtxt(mag_data_path)
                    
                    # Handle single-line files
                    if data.ndim == 1:
                        data = data.reshape(1, -1)
                    
                    # Extract temperature and magnetisation columns
                    temperature, magnetisation = data[:, 0], data[:, 1]
                    
                    # Sort by temperature
                    sorted_indices = np.argsort(temperature)
                    temperature_sorted = temperature[sorted_indices]
                    magnetisation_sorted = magnetisation[sorted_indices]
                    
                    # Store in dictionary
                    data_dict[N_value] = (temperature_sorted, magnetisation_sorted)
                    print(f"Loaded data for N={N_value}")
                else:
                    print(f"Warning: No 'mag_data' file found in {item_path}")
            except ValueError:
                print(f"Warning: Could not parse N value from folder name '{item}'")
    
    return data_dict

def plot_magnetisation(data_dict):
    if not data_dict:
        print("No data to plot.")
        return
        
    plt.figure(figsize=(10, 7))
    # Define marker styles, colors, and line styles for different plots
    markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'P', 'h']
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']
    line_styles = ['-', '--', '-.', ':']
    
    
    # Plot each dataset
    for i, N in enumerate(N_values):
        temperature, magnetisation = data_dict[N]
        plt.plot(temperature, magnetisation, marker=markers[i], linestyle=line_styles[i], c=colors[i], label=f"N={N}")
    
    plt.xlabel("Temperature", fontsize=12)
    plt.ylabel("Avg Magnetisation per Particle", fontsize=12)
    plt.title("Temperature vs. Magnetisation for Different N Values", fontsize=14)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("magnetisation_plot.png", dpi=300)  # Save the figure
    plt.show()

if __name__ == "__main__":
    data_dict = read_magnetisation_data()
    plot_magnetisation(data_dict)
