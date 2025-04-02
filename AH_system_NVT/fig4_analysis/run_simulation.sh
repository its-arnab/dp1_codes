#!/bin/bash

SOURCE_FILE="ah_nvt_cell.c"
EXECUTABLE_BASE="ah_nvt_cell"

N_VALUES=(200 2000)
K_VALUES=(0.1 1.0)

# Generate log-spaced temperatures using awk
TEMPERATURES=($(awk 'BEGIN {
    n = 20; min_exp = -2; max_exp = 0.8;
    for (i = 0; i < n; i++) {
        val = 10^(min_exp + i * (max_exp - min_exp) / (n - 1));
        printf("%.5f ", val);
    }
}'))

### Prepare directories ###
mkdir -p log_files avg_data Excess_data
rm -rf log_files/* avg_data/* Excess_data/*

for N in "${N_VALUES[@]}"; do
    for K in "${K_VALUES[@]}"; do
        K_str=$(printf "%.1f" "$K" | tr '.' '_')  # 0.1 → 0_1
        for T in "${TEMPERATURES[@]}"; do
            T_str=$(printf "%.5f" "$T" | tr '.' '_')  # 0.01 → 0_01

            # Create unique C file and executable for each parameter set
            C_FILE="ah_nvt_cell_cp_N${N}_K${K_str}_T${T_str}.c"
            EXECUTABLE="${EXECUTABLE_BASE}_N${N}_K${K_str}_T${T_str}.out"
            
            # Copy and modify the C file
            cp "$SOURCE_FILE" "$C_FILE"
            sed -i "15s/^#\s*define\s*N\s*.*/#define N $N/" "$C_FILE"
            sed -i "24s/^const double K = .*/const double K = $K;/" "$C_FILE"
            sed -i "21s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"

            # Compile, run, and clean up after completion
            gcc -o "$EXECUTABLE" "$C_FILE" -lm
            (nohup "./$EXECUTABLE" > "log_files/output_N_${N}_K_${K_str}_T_${T_str}.log" 2>&1 \
              && rm "$EXECUTABLE" "$C_FILE") &
        done
    done
done