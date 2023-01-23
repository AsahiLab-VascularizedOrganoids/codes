#!/bin/bash

####### Params ########################
target_dir="$HOME/resources/hETV2/"
sample_names=(
    "hETV2_hCO" "hETV2_vhCO"
)
reference_genome="${HOME}/ref/refdata-gex-GRCh38-2020-A"
#######################################

mkdir -p "${target_dir}/counted"
cd "${target_dir}" || exit
bamtofastq --nthreads=4 "raw" "fq"

for sample_name in "${sample_names[@]}"; do
    cellranger count --id="${sample_name}" --sample="${sample_name}" --fastqs="$target_dir/fq" \
        --transcriptome="$reference_genome" --no-bam
done
