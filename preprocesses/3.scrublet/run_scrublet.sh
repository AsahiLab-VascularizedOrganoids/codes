#!/bin/bash
sample_names=(
    "Fetal" "hETV2_hCO" "hETV2_vhCO"
    "HUVEC_hCO_d65" "HUVEC_hCO_d100"
    "HUVEC_vhCO_d65" "HUVEC_vhCO_d100"
)
root_dir="$HOME/resources/"

for sample_name in "${sample_names[@]}"; do
    input="${root_dir}/${sample_name}/counted/counts.tsv"
    output="${root_dir}/${sample_name}/counted/filtered_counts.tsv"
    python "./functions/scrublet.py" "$sample_name" "$input" "$output"
done
