#!/bin/bash
#PBS -N KA_T_$TEMP
#PBS -A project_code
#PBS -l walltime=72:00:00
#PBS -q cpuq
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=32

cd $PBS_O_WORKDIR

# Define temperature and number of ensemble runs
TEMP=$1
N_ENSEMBLE=32


# Convert TEMP to delay seconds: e.g., (TEMP - 0.50) * 200
# To avoid floating point issues in bash, use bc

DELAY=$(echo "($TEMP - 0.50) * 200" | bc)
echo "Sleeping for $DELAY seconds to avoid conflicts..."
sleep $DELAY


# Modify the global.h to set the correct temperature
sed -i "72s/^.*$/#define     T_res       $TEMP/" include/global.h


make || { echo "Build failed"; exit 1; }

# Log directory
mkdir -p log_files/T_$TEMP/

# Run all ensembles
for i in $(seq 0 $((N_ENSEMBLE - 1))); do
    run_id=$(printf "run_%02d" $i)
    echo ">> Starting $run_id for T=$TEMP..."
    nohup ./a_exec new "$run_id" > "log_files/T_$TEMP/log_${run_id}.txt" 2>&1 &
done

wait
