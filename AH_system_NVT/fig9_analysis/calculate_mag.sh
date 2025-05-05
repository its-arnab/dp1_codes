#!/bin/bash

# Directory and output settings
N="10000"
DIR="magnetisation_data/N_${N}"
ENSEMBLE_PREFIX="${DIR}/run_"
OUTPUT_FILE="${DIR}/mag_data"


# Temporary file to store concatenated and sorted data
TMP_SORTED="temp_sorted_data.txt"
> "$TMP_SORTED"  # Clear if exists

# Concatenate and sort data from all run_XX files
for i in $(seq -w 0 15); do
    FILE="${ENSEMBLE_PREFIX}${i}"
    if [ -f "$FILE" ]; then
        awk '{printf "%.5f %.8f\n", $1, $2}' "$FILE" >> "$TMP_SORTED"
    else
        echo "Warning: File $FILE not found, skipping..."
    fi
done

# Sort the data numerically by temperature (1st column)
sort -n "$TMP_SORTED" > "${TMP_SORTED}.sorted"

# Use awk to compute average magnetization grouped by temperature
awk '
{
    temp = $1
    mag = $2
    sum[temp] += mag
    count[temp]++
}
END {
    for (t in sum) {
        printf "%.5f %.8f\n", t, sum[t] / count[t]
    }
}
' "${TMP_SORTED}.sorted" | sort -n > "$OUTPUT_FILE"

# Clean up
rm -f "$TMP_SORTED" "${TMP_SORTED}.sorted"

echo "Averaged data saved to: $OUTPUT_FILE"

