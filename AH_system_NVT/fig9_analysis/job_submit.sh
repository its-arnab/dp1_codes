#!/bin/bash

# ------------------------- Temperature Ranges -------------------------
# TEMP1=($(awk 'BEGIN {for(i=0; i<5; i++) printf "%.5f ", 0.01 + i * (0.30 - 0.01) / 4}'))
# TEMP2=($(awk 'BEGIN {for(i=0; i<10; i++) printf "%.5f ", 0.35 + i * (0.50 - 0.35) / 9}'))
# TEMP3=($(awk 'BEGIN {for(i=0; i<5; i++) printf "%.5f ", 0.60 + i * (1.00 - 0.60) / 4}'))

# TEMPERATURES=("${TEMP1[@]}" "${TEMP2[@]}" "${TEMP3[@]}")

TEMPERATURES=(0.01)

# ------------------------- Template Check -------------------------
TEMPLATE_RUN="run.sh"
if [ ! -f "$TEMPLATE_RUN" ]; then
    echo "Error: Template $TEMPLATE_RUN not found."
    exit 1
fi

chmod +x "$TEMPLATE_RUN"
mkdir -p log_files
rm -f log_files/*

# ------------------------- Submission -------------------------
for T in "${TEMPERATURES[@]}"; do
    T_clean=$(echo "$T" | sed 's/\.//')
    RUN_SCRIPT="run_T_${T_clean}.sh"

    # Create customized run script
    sed "s/__TEMP__/$T/g" "$TEMPLATE_RUN" > "$RUN_SCRIPT"
    chmod +x "$RUN_SCRIPT"

    echo "Submitting job for T = $T"
    nohup bash "$RUN_SCRIPT" > "log_files/log_master_T_${T}.log" 2>&1 &
    wait
    rm "$RUN_SCRIPT"
    echo "Job for T = $T completed and cleaned up."
done

chmod -x "$TEMPLATE_RUN"
echo "All temperature jobs submitted."
