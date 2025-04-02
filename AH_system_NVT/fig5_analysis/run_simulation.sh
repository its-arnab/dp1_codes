#!/bin/bash

SOURCE_FILE="ah_nvt_cell.c"
cp "$SOURCE_FILE" ah_nvt_cell_cp.c      # make a copy
C_FILE="ah_nvt_cell_cp.c"
EXECUTABLE="ah_nvt_cell.out"

T_VALUES=(0.3 0.4 1.0)
K_VALUES=(0.1 1.0)

### Clean up every folder ###
mkdir -p log_files movie_data
rm -rf log_files/* movie_data/*
### ===================== ###

for T in "${T_VALUES[@]}"; do
    T_str=$(printf "%.5f" "$T" | tr '.' '_')  # 0.01 → 0_01
    for K in "${K_VALUES[@]}"; do
        K_str=$(printf "%.1f" "$K" | tr '.' '_')  # 0.1 → 0_1

        # Create unique C file and executable for each parameter set
        C_FILE="ah_nvt_cell_cp_K${K_str}_T${T_str}.c"
        EXECUTABLE="${EXECUTABLE_BASE}_K${K_str}_T${T_str}.out"

        # Copy and modify the C file
        cp "$SOURCE_FILE" "$C_FILE"

        # Modify source code
        sed -i "21s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"
        sed -i "24s/^const double K = .*/const double K = $K;/" "$C_FILE"

        # Compile, run, and clean up after completion
        gcc -o "$EXECUTABLE" "$C_FILE" -lm
        (nohup "./$EXECUTABLE" > "log_files/output_N_${N}_K_${K_str}_T_${T_str}.log" 2>&1 \
            && rm "$EXECUTABLE" "$C_FILE") &
    done
done
