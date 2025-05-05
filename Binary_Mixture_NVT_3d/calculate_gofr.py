import os
import numpy as np
from glob import glob

# ===== USER PARAMETERS ===== #
base_path   = "../data_files"
temperature = "T_1.10000"
runs        = [f"run_{i:02d}" for i in range(32)]
blocksize   = 80000

# — Physical parameters —
N   = 1000         # number of particles of type‑0
rho = 1.20         # number density
L   = (N / rho)**(1/3)    # box length
r_max  = L / 2
n_bins  = 500             # you may adjust for resolution

dr = r_max / n_bins

# radii for shell centers
edges = np.linspace(0, r_max, n_bins+1)
r_centers = 0.5 * (edges[:-1] + edges[1:])
shell_vol = 4 * np.pi * r_centers**2 * dr
nid = shell_vol * rho      # number of particles in the shell(ideal gas)


def read_positions(filename):
    data = np.loadtxt(filename)
    mask = (data[:,0] == 0)
    pos  = data[mask, 1:4]    # shape (N,3)
    pos %= L   # fold into [0,L)
    return pos

def compute_gr(positions):
    """
    Compute g(r) for a single snapshot of unwrapped positions:
    - positions: (N,3) array, already folded into [0,L)
    Returns (r_centers, g_r)
    """
    # pair histogram
    hist = np.zeros(n_bins, dtype=int)
    Np = positions.shape[0]

    for i in range(Np-1):
        delta = positions[i+1:] - positions[i]
        # minimum-image PBC
        delta -= L * np.round(delta / L)
        dist2  = np.sum(delta**2, axis=1)
        dist   = np.sqrt(dist2)
        bins   = np.floor(dist / dr).astype(int)
        valid = (bins >= 0) & (bins < n_bins)
        np.add.at(hist, bins[valid], 2)

    # normalization by N^2
    norm = Np * nid
    g_r = hist / norm
    return  g_r

def compute_gr_for_run(run_path):
    # collect and sort available step numbers
    files     = sorted(glob(os.path.join(run_path, "movie_data", "step_*")))
    step_nums = sorted(int(os.path.basename(f).split("_")[-1]) for f in files)
    if not step_nums:
        print(f"  [run={os.path.basename(run_path)}] no data found")
        return None

    # last block start
    block_starts = [s for s in range(1, step_nums[-1]+1, blocksize)]
    last_start   = block_starts[-1]

    gr_accum = []
    for n in range(15):
        step = (last_start - 1) + (1 << n)
        if step not in step_nums:
            print(f"  [run={os.path.basename(run_path)}] missing step {step}, skipping")
            continue
        fname = os.path.join(run_path, "movie_data", f"step_{step:07d}")
        pos   = read_positions(fname)
        g  = compute_gr(pos)
        gr_accum.append(g)

    if not gr_accum:
        return None

    # average over the four frames
    avg_gr = np.mean(np.vstack(gr_accum), axis=0)
    return avg_gr  # return avg_gr

def process_temperature_gr(temp_path):
    all_gr = []
    for rname in runs:
        res = compute_gr_for_run(os.path.join(temp_path, rname))
        if res is not None:
            all_gr.append(res)

    if not all_gr:
        print("No g(r) data found across all runs.")
        return

    ensemble_gr = np.mean(np.vstack(all_gr), axis=0)

    out = os.path.join(temp_path, f"gr_{os.path.basename(temp_path)}.dat")
    header = f"# r  g(r)   # L={L:.5f}, n_bins={n_bins}"
    np.savetxt(out, np.column_stack((r_centers, ensemble_gr)),
               header=header, fmt="%.8e")
    print(f"Saved ensemble g(r): {out}")

if __name__ == "__main__":
    temp_path = os.path.join(base_path, temperature)
    process_temperature_gr(temp_path)
