#!/bin/bash
#PBS -N ah_N2000
#PBS -A project_code
#PBS -l walltime=4:00:00
#PBS -q cpuq
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=4

cd $PBS_O_WORKDIR

# Create log directory
mkdir -p log_files


SOURCE_FILE="ah_nvt_cell.c"
EXECUTABLE_BASE="ah_nvt_cell"
T_VALUES=(0.02 0.2 0.4 0.5)

for T in "${T_VALUES[@]}"; do
    T_str=$(printf "%.5f" "$T" | tr '.' '_')  # 0.01 → 0_01

    # Create unique C file and executable for each parameter set
    C_FILE="ah_nvt_cell_cp_T${T_str}.c"
    EXECUTABLE="${EXECUTABLE_BASE}_T${T_str}.out"

    # Copy and modify the C file
    cp "$SOURCE_FILE" "$C_FILE"

    # Modify source code (ensure line 21 is correct)
    sed -i "21s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"

    # Compile, run, and clean up after completion
    gcc -o "$EXECUTABLE" "$C_FILE" -lm
    (nohup "./$EXECUTABLE" > "log_files/output_T_${T_str}.log" 2>&1 \
        && rm "$EXECUTABLE" "$C_FILE") &
done

wait