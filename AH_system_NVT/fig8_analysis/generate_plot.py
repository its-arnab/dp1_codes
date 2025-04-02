import os
import re
import numpy as np
import matplotlib.pyplot as plt

def extract_parameters(filename):
    """Extracts T from the filename."""
    match = re.search(r'T_([0-9]+\.[0-9]+)', filename)
    if match:
        T = float(match.group(1))
        return T
    return None

def read_velocity_data(filepath):
    """Reads vx data from file."""
    try:
        data = np.loadtxt(filepath, usecols=(0, 1))  # Assuming space-separated values
        vx = data[:, 0]
        return vx
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return np.array([])

def theoretical_gaussian(vx, T):
    """Computes theoretical Gaussian distribution P_th(v_x)."""
    return np.sqrt(1 / (2 * np.pi * T)) * np.exp(-vx**2 / (2 * T))

def compute_ensemble_histograms(velocity_path, bins=50):
    """Computes ensemble-averaged histograms for different temperatures."""
    temp_data = {}
    
    print("Checking available runs...")
    
    for run in range(21):  # Loop over runs run_00 to run_20
        run_folder = os.path.join(velocity_path, f'run_{run:02d}')
        if not os.path.isdir(run_folder):
            print(f"Run folder missing: {run_folder}")
            continue  # Skip if folder does not exist
        
        print(f"Processing {run_folder}...")
        
        for filename in os.listdir(run_folder):
            T = extract_parameters(filename)
            if T is None:
                print(f"Skipping file (T not found): {filename}")
                continue  # Skip if filename format is incorrect
            
            filepath = os.path.join(run_folder, filename)
            vx = read_velocity_data(filepath)
            
            if vx.size == 0:
                print(f"Skipping empty file: {filepath}")
                continue
            
            hist, bin_edges = np.histogram(vx, bins=bins, density=True)
            
            if T not in temp_data:
                temp_data[T] = [hist, bin_edges]
            else:
                temp_data[T][0] += hist  # Sum histograms
    
    print("\nAveraging histograms...")
    for T in temp_data:
        temp_data[T][0] /= 21  # Normalize by number of runs
        print(f"T = {T}: Histogram averaged over 21 runs.")
    
    return temp_data

def plot_histograms(temp_data, selected_temps):
    """Plots ensemble-averaged histograms for selected temperatures."""
    plt.figure(figsize=(8, 6))
    
    for T in selected_temps:
        if T not in temp_data:
            print(f"Warning: T = {T} not found in computed data!")
            continue
        
        hist, bin_edges = temp_data[T]
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute middle points
        plt.scatter(bin_centers, hist, label=f'simu T = {T}')
        plt.yscale('log')
        
        # Plot theoretical Gaussian
        vx_range = np.linspace(bin_centers.min(), bin_centers.max(), 200)
        P_th = theoretical_gaussian(vx_range, T)
        plt.plot(vx_range, P_th, linestyle='-', label=f'Th T = {T}')
    
    plt.xlabel('$v_x$', fontsize=14)
    plt.ylabel('$P(v_x)$', fontsize=14)
    plt.legend()
    plt.title('Ensemble Averaged $P(v_x)$ vs $v_x$', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

def main():
    velocity_path = "Vcm_data/"
    temp_data = compute_ensemble_histograms(velocity_path)
    
    print("\nExtracted temperatures:", sorted(temp_data.keys()))
    selected_temps = sorted(temp_data.keys())[:3]  # Pick three temperatures
    
    if len(selected_temps) < 3:
        print("Warning: Less than 3 temperatures found!")
    
    plot_histograms(temp_data, selected_temps)

if __name__ == "__main__":
    main()
