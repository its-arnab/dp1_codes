import os
import numpy as np
import matplotlib.pyplot as plt

def read_magnetisation_data(folder="magnetisation_data"):
    data_dict = {}
    
    for filename in os.listdir(folder):
        if filename.startswith("N_") and filename[2:].isdigit():
            N_value = int(filename[2:])  # Extracting N from filename
            filepath = os.path.join(folder, filename)
            
            # Read file content
            data = np.loadtxt(filepath)
            if data.ndim == 1:
                data = data.reshape(1, -1)  # Handle single-line files
            
            temperature, magnetisation = data[:, 0], data[:, 1]
            
            # Sort by temperature
            sorted_indices = np.argsort(temperature)
            temperature_sorted = temperature[sorted_indices]
            magnetisation_sorted = magnetisation[sorted_indices]
            
            data_dict[N_value] = (temperature_sorted, magnetisation_sorted)
    
    return data_dict

def plot_magnetisation(data_dict):
    plt.figure(figsize=(8, 6))
    
    for N, (temperature, magnetisation) in sorted(data_dict.items()):
        plt.plot(temperature, magnetisation, marker='o', linestyle='-', label=f"N={N}")
    
    plt.xlabel("Temperature")
    plt.ylabel("Avg Magnetisation per Particle")
    plt.title("Temperature vs. Magnetisation for Different N Values")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    data_dict = read_magnetisation_data()
    plot_magnetisation(data_dict)
