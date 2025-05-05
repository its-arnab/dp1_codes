import os
import numpy as np
from glob import glob

# ===== USER PARAMETERS =====
dt = 0.005
base_path = "../data_files"
temperature = "T_1.10000"
runs = [f"run_{i:02d}" for i in range(32)]
blocksize = 80000
k_mag = 7.25

# ===== FUNCTION TO READ POSITION FILE =====
def read_positions(filename):
    data = np.loadtxt(filename)
    type_0 = data[:, 0] == 0
    return data[type_0][:, 1:]  # Only x, y, z of type-0

# ===== Fs(k,t) CALCULATION FUNCTION FOR A SINGLE RUN =====
def compute_fskt_for_run(run_path):
    files = sorted(glob(os.path.join(run_path, "movie_data", "step_*")))
    files.sort(key=lambda x: int(os.path.basename(x).split("_")[-1]))

    if not files:
        return None

    step_nums = [int(os.path.basename(f).split("_")[-1]) for f in files]
    step_to_index = {step: idx for idx, step in enumerate(step_nums)}

    traj = [read_positions(f) for f in files]
    traj = np.array(traj)  # shape: (n_frames, N0, 3)

    fskt_data = []

    # ===== Intra-block Analysis =====
    for block_start in range(1, step_nums[-1] + 1, blocksize):
        origin_step = block_start

        if origin_step not in step_to_index:
            print(f"Block start {block_start} not in step_to_index, skipping.")
            continue

        origin_idx = step_to_index[origin_step]
        r0 = traj[origin_idx]

        for n in range(1, 17):  # 2^1 to 2^16
            lag_step = (origin_step - 1) + (1 << n)

            if lag_step > step_nums[-1] or lag_step not in step_to_index:
                continue

            lag_idx = step_to_index[lag_step]
            r_t = traj[lag_idx]
            dx = r_t[:, 0] - r0[:, 0]  # Only x-component
            cos_kdx = np.cos(k_mag * dx)
            fskt = np.mean(cos_kdx) 
            fskt_data.append((lag_step - origin_step, fskt))

    # ===== Inter-block Analysis =====
    block_starts = [s for s in range(1, step_nums[-1] + 1, blocksize) if s in step_to_index]

    for i, origin_step in enumerate(block_starts[:-1]):
        origin_idx = step_to_index[origin_step]
        r0 = traj[origin_idx]

        for target_start in block_starts[i + 1:]:
            for n in range(1, 17):
                target_step = (target_start - 1) + (1 << n)
                if target_step in step_to_index:
                    target_idx = step_to_index[target_step]
                    r_t = traj[target_idx]
                    dx = r_t[:, 0] - r0[:, 0]  # Only x-component
                    cos_kdx = np.cos(k_mag * dx)
                    fskt = np.mean(cos_kdx)
                    fskt_data.append((target_step - origin_step, fskt))
                else:
                    continue

    if not fskt_data:
        return None

    # Group by lag time and average
    fskt_data = np.array(fskt_data)
    lag_times = fskt_data[:, 0]
    f_vals = fskt_data[:, 1]
    unique_lags = np.unique(lag_times)
    avg_fskt = np.array([f_vals[lag_times == t].mean() for t in unique_lags])
    unique_lags *= dt

    return unique_lags, avg_fskt

# ===== MAIN SEQUENTIAL EXECUTION FOR ALL RUNS =====
def process_temperature(temp_path):
    results = []
    for r in runs:
        run_path = os.path.join(temp_path, r)
        result = compute_fskt_for_run(run_path)
        if result is not None:
            results.append(result)

    if not results:
        print("No valid runs found.")
        return

    time_ref = results[0][0]
    fskt_vals = np.array([res[1] for res in results])
    avg_fskt = np.mean(fskt_vals, axis=0)

    output_file = os.path.join(temp_path, f"fskt_k{int(k_mag*100):05d}_{os.path.basename(temp_path)}.dat")
    np.savetxt(output_file, np.column_stack((time_ref, avg_fskt)), header="time Fs(k,t)", fmt="%.8e")
    print(f"Saved Fs(k,t): {output_file}")

# ===== RUN THE SCRIPT FOR ONE TEMPERATURE =====
if __name__ == "__main__":
    temp_path = os.path.join(base_path, temperature)
    process_temperature(temp_path)
