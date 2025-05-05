#!/bin/bash

N_ENSEMBLE=32  # You can change this as needed
LOG_DIR="log_files/T_1.10"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Build once
make || { echo "Build failed"; exit 1; }

# Run each ensemble in background using taskset + nohup
for i in $(seq 16 $((N_ENSEMBLE - 1))); do
    run_id=$(printf "run_%02d" $i)

    echo ">> Starting $run_id on core $core..."
    nohup ./a_exec new "$run_id" > "${LOG_DIR}/log_${run_id}.txt" 2>&1 &
done

echo "Waiting for all background processes to finish..."
wait  # Wait for all jobs to complete

# Clean up executables
echo "Cleaning up build files..."
make clean

echo "All ensemble runs completed. Logs are in '${LOG_DIR}' and build files cleaned."
