#!/bin/bash

# ------------------------- Configuration -------------------------
SOURCE_FILE="ah_nvt_cell.c"
EXECUTABLE_BASE="ah_nvt_cell"
T=__TEMP__      # <<< THIS LINE will be replaced in master script

# ------------------------- Naming Conventions -------------------------
T_str=$(printf "%.5f" "$T" | tr '.' '_')
C_FILE="ah_nvt_cell_cp_T${T_str}.c"
EXECUTABLE="${EXECUTABLE_BASE}_T${T_str}.out"

# ------------------------- Ensemble Generation -------------------------
ENSEMBLES=()
for i in $(seq -w 0 16); do
    ENSEMBLES+=("run_$i")
done

# ------------------------- File Preparation -------------------------
cp "$SOURCE_FILE" "$C_FILE"

sed -i "22s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"

gcc -o "$EXECUTABLE" "$C_FILE" -lm -O2

mkdir -p log_files

# ------------------------- Run Simulations -------------------------
for RUN_ID in "${ENSEMBLES[@]}"; do
    (nohup "./$EXECUTABLE" "$RUN_ID" > "log_files/output_${RUN_ID}.log" 2>&1) &
done

wait

# ------------------------- Clean-up -------------------------
rm "$EXECUTABLE" "$C_FILE"
