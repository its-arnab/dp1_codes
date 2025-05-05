import os
import numpy as np
from glob import glob
import _multiprocessing as mp

# ===== USER PARAMETERS =====
dt = 0.005
base_path = "../data_files"
temperature = "T_0.50000"
runs = [f"run_{i:02d}" for i in range(16)]
blocksize = 80000
a = 0.3
a2 = a * a

# ===== FUNCTION TO READ POSITION FILE =====
def read_positions(filename):
    data = np.loadtxt(filename)
    type_0 = data[:, 0] == 0
    return data[type_0][:, 1:]  # Only x, y, z of type-0

# ===== Q(t) CALCULATION FUNCTION FOR A SINGLE RUN =====
def compute_qoft_for_run(run_path):
    files = sorted(glob(os.path.join(run_path, "movie_data", "step_*")))
    files.sort(key=lambda x: int(os.path.basename(x).split("_")[-1]))

    if not files:
        return None

    step_nums = [int(os.path.basename(f).split("_")[-1]) for f in files]
    step_to_index = {step: idx for idx, step in enumerate(step_nums)}

    # print(step_to_index)

    traj = [read_positions(f) for f in files]
    traj = np.array(traj)  # shape: (n_frames, N0, 3)

    qoft_data = []

    # ===== Intra-block Analysis =====
    for block_start in range(1, step_nums[-1] + 1, blocksize):
        origin_step = block_start

        if origin_step not in step_to_index:
            print(f"Block start {block_start} not in step_to_index, skipping.")   # report problem
            continue

        origin_idx = step_to_index[origin_step]
        for n in range(1, 17):  # 2^0 to 2^16 = 65536
            lag_step = (origin_step - 1) + (1 << n)
    
            if lag_step > step_nums[-1] or lag_step not in step_to_index:
                print(f"Processing intra-block comparison for block starting at step {block_start}, lag step {lag_step}")
                continue

            lag_idx = step_to_index[lag_step]
            disp = traj[lag_idx] - traj[origin_idx]
            squared = np.sum(disp**2, axis=1)
            qoft = np.sum(squared < a2) / len(squared)
            qoft_data.append((lag_step - origin_step, qoft))

    # ===== Inter-block Analysis (Pairwise Block Comparisons) =====
    block_starts = [s for s in range(1, step_nums[-1] + 1, blocksize) if s in step_to_index]

    # print(f"Block starts: {block_starts}")
    
    for i, origin_step in enumerate(block_starts[:-1]):
        # print(f"Processing inter-block comparison for block starting at step {i, origin_step}")
        origin_idx = step_to_index[origin_step]
        for target_start in block_starts[i + 1:]:
            for n in range(1, 17):  # lag steps: 2^0 to 2^16
                target_step = (target_start - 1) + (1 << n)
                if target_step in step_to_index:
                    target_idx = step_to_index[target_step]
                    disp = traj[target_idx] - traj[origin_idx]
                    squared = np.sum(disp**2, axis=1)
                    qoft = np.sum(squared < a2) / len(squared)
                    qoft_data.append((target_step - origin_step, qoft))
                else:
                    print(f"Target step {target_step} not in step_to_index, skipping.")
                    continue


    if not qoft_data:
        return None

    # Group by lag time and average
    qoft_data = np.array(qoft_data)
    lag_times = qoft_data[:, 0]
    q_vals = qoft_data[:, 1]
    unique_lags = np.unique(lag_times)
    avg_qoft = np.array([q_vals[lag_times == t].mean() for t in unique_lags])
    unique_lags *= dt

    return unique_lags, avg_qoft


# ===== MAIN SEQUENTIAL EXECUTION FOR ALL RUNS =====
def process_temperature(temp_path):

    results = []
    for r in runs:
        run_path = os.path.join(temp_path, r)
        result = compute_qoft_for_run(run_path)
        if result is not None:
            results.append(result)

    if not results:
        print("No valid runs found.")
        return

    time_ref = results[0][0]
    qoft_vals = np.array([res[1] for res in results])
    avg_qoft = np.mean(qoft_vals, axis=0)

    output_file = os.path.join(temp_path, f"qoft_{os.path.basename(temp_path)}.dat")
    np.savetxt(output_file, np.column_stack((time_ref, avg_qoft)), header="time Q(t)", fmt="%.8e")
    print(f"Saved Q(t): {output_file}")

# ===== RUN THE SCRIPT FOR ONE TEMPERATURE =====
if __name__ == "__main__":
    temp_path = os.path.join(base_path, temperature)
    process_temperature(temp_path)


















'''
import os
import numpy as np
from glob import glob
import multiprocessing as mp

# ===== USER PARAMETERS =====
dt = 0.005
base_path = "../data_files"
temperature = "T_0.50000"
runs = [f"run_{i:02d}" for i in range(8)]
save_folder = "msd_results"
os.makedirs(save_folder, exist_ok=True)

# ===== FUNCTION TO READ POSITION FILE =====
def read_positions(filename):
    data = np.loadtxt(filename)
    type_0 = data[:, 0] == 0
    return data[type_0][:, 1:]  # Only x, y, z of type-0

# ===== MSD CALCULATION FUNCTION FOR A SINGLE RUN =====
def compute_msd_for_run(run_path):
    files = sorted(glob(os.path.join(run_path, "movie_data", "step_*")))
    files.sort(key=lambda x: int(os.path.basename(x).split("_")[-1]))

    n_frames = len(files)
    if n_frames == 0:
        return None

    traj = [read_positions(f) for f in files]
    traj = np.array(traj)  # shape: (n_frames, N0, 3)
    N0 = traj.shape[1]

    max_lag = int(np.log2(n_frames))
    msd_list = []
    time_list = []

    for k in range(1, max_lag):
        lag = 2**k
        displacements = traj[lag:] - traj[:-lag]
        squared = np.sum(displacements**2, axis=2)
        msd = np.mean(squared)
        msd_list.append(msd)
        time_list.append(lag * dt)

    return np.array(time_list), np.array(msd_list)

# ===== MAIN PARALLEL EXECUTION FOR ALL RUNS =====
def process_temperature(temp_path):
    pool = mp.Pool(processes=mp.cpu_count())
    run_paths = [os.path.join(temp_path, r) for r in runs]
    results = pool.map(compute_msd_for_run, run_paths)
    pool.close()
    pool.join()

    # Remove None entries (if any run had no data)
    results = [r for r in results if r is not None]
    if not results:
        print("No valid runs found.")
        return

    time_ref = results[0][0]
    msds = np.array([res[1] for res in results])
    avg_msd = np.mean(msds, axis=0)

    output_file = os.path.join(save_folder, f"msd_{os.path.basename(temp_path)}.dat")
    np.savetxt(output_file, np.column_stack((time_ref, avg_msd)), header="time MSD", fmt="%.8e")
    print(f"Saved MSD: {output_file}")

# ===== RUN THE SCRIPT FOR ONE TEMPERATURE =====
if __name__ == "__main__":
    temp_path = os.path.join(base_path, temperature)
    process_temperature(temp_path)

    

'''
