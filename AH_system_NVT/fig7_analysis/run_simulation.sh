SOURCE_FILE="ah_nvt_cell.c"
T_VALUES=(1.0 2.0 3.0)

for T in "${T_VALUES[@]}"; do
    T_str=$(printf "%.5f" "$T" | tr '.' '_')  # 0.01 → 0_01
    
    for run in {0..20}; do
        # Create unique C file and executable for each parameter set
        C_FILE="ah_nvt_cell_cp_T${T_str}_run${run}.c"
        EXECUTABLE="${EXECUTABLE_BASE}_T${T_str}_run${run}.out"

        # Copy and modify the C file
        cp "$SOURCE_FILE" "$C_FILE"

        # Modify source code
        sed -i "21s/^const double T_res = .*/const double T_res = $T;/" "$C_FILE"
        sed -i "86s/^int run = .*/int run = $run;/" "$C_FILE"

        # Compile, run, and clean up after completion
        gcc -o "$EXECUTABLE" "$C_FILE" -lm
        (nohup "./$EXECUTABLE" > "log_files/output_T_${T_str}_run_${run}.log" 2>&1 \
            && rm "$EXECUTABLE" "$C_FILE") &
    done

done
