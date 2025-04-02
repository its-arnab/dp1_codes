#!/bin/bash

SOURCE_FILE="ah_nvt_cell.c"
EXECUTABLE_BASE="ah_nvt_cell"

N_VALUES=(200)

# Generate first range: 15 points from 0.05 to 0.5
TEMP1=($(awk 'BEGIN {for(i=0; i<15; i++) printf "%.5f ", 0.05 + i * (0.5 - 0.05) / 14}'))

# Generate second range: 5 points from 0.6 to 1.0
TEMP2=($(awk 'BEGIN {for(i=0; i<5; i++) printf "%.5f ", 0.6 + i * (1.0 - 0.6) / 4}'))

# Combine both arrays
TEMPERATURES=("${TEMP1[@]}" "${TEMP2[@]}")

### ==== Prepare directories ==== ###
mkdir -p log_files magnetisation_data final_config
rm -rf log_files/* magnetisation_data/* final_config/*
### ============================= ###

for N in "${N_VALUES[@]}"; do
    for T in "${TEMPERATURES[@]}"; do
        T_str=$(printf "%.5f" "$T" | tr '.' '_')  # 0.01 → 0_01

        # Create unique C file and executable for each parameter set
        C_FILE="ah_nvt_cell_cp_N${N}_T${T_str}.c"
        EXECUTABLE="${EXECUTABLE_BASE}_N${N}_T${T_str}.out"
        
        # Copy and modify the C file
        cp "$SOURCE_FILE" "$C_FILE"
        sed -i "15s/^#\s*define\s*N\s*.*/#define N $N/" "$C_FILE"
        sed -i "21s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"

        # Compile, run, and clean up after completion
        gcc -o "$EXECUTABLE" "$C_FILE" -lm
        (nohup "./$EXECUTABLE" > "log_files/output_N_${N}_T_${T_str}.log" 2>&1 \
            && rm "$EXECUTABLE" "$C_FILE") &
    done
done
