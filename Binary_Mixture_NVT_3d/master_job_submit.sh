#!/bin/bash

# List of temperatures
temps=(0.50 0.60 0.70 0.80 0.90 1.00)

for T in "${temps[@]}"; do
    echo "Submitting job for T = $T"
    qsub -v TEMP=$T run_k4.sh
done
