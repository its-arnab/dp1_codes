#!/bin/bash

echo "-------------------------------------------"

output_file="../data_files/thermo_avg"
> "$output_file"  # clear previous contents

for T in 0.44000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000 1.10000; do
    temp_path="../data_files/T_${T}"

    # Temporary file to store per-run averages
    tmp_avg_file="${temp_path}/tmp_avg_data"
    > "$tmp_avg_file"

    # Loop over ensemble runs
    for run in $(seq -w 0 31); do
        file="${temp_path}/run_${run}/T_P_KE_PE_E_data"
        if [[ -f "$file" ]]; then
            tail -n +100001 "$file" | \
            awk '{
                for(i=1;i<=5;i++) sum[i]+=$i;
                count++
            }
            END {
                if (count > 0) {
                    for(i=1;i<=5;i++) printf "%.6f ", sum[i]/count;
                    print ""
                }
            }' >> "$tmp_avg_file"
        fi
    done

    # Average over all runs
    if [[ -s "$tmp_avg_file" ]]; then
        awk -v Tval="$T" '{
            for(i=1;i<=NF;i++) sum[i]+=$i;
            count++
        }
        END {
            printf "%.5f ", Tval;
            for(i=1;i<=NF;i++) printf "%.6f ", sum[i]/count;
            print ""
        }' "$tmp_avg_file" >> "$output_file"
        rm "$tmp_avg_file"
    else
        echo "Warning: No valid data found for T = $T"
    fi
done

echo "-------------------------------------------"
