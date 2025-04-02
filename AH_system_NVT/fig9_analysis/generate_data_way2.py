import os
import numpy as np

def compute_magnetisation(folder="final_config", output_folder="magnetisation_data_way2"):
    os.makedirs(output_folder, exist_ok=True)  # Ensure output directory exists
    
    for filename in os.listdir(folder):
        if filename.startswith("N_") and "_T_" in filename:
            try:
                parts = filename.split("_T_")
                N_value = int(parts[0][2:])  # Extract N
                T_value = float(parts[1])    # Extract T
                filepath = os.path.join(folder, filename)
                
                # Read file content
                data = np.loadtxt(filepath)                
                theta = data[:, 2]
                mx = np.sum(np.cos(theta)) / N_value
                my = np.sum(np.sin(theta)) / N_value
                m = np.sqrt(mx**2 + my**2)
                
                # Prepare output file
                output_file = os.path.join(output_folder, f"N_{N_value:05d}")
                with open(output_file, "a") as f:
                    f.write(f"{T_value:.5f} {mx:.6f} {my:.6f} {m:.6f}\n")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                continue

if __name__ == "__main__":
    compute_magnetisation()
