#!/bin/bash

####### Params ########################
target_dir="$HOME/resources/HUVEC/"
sample_names=(
    "HUVEChCOd65" "HUVECvhCOd65"
    "HUVEChCOd100" "HUVECvhCOd100"
)
reference_genome="${HOME}/ref/refdata-gex-GRCh38-2020-A"
#######################################

mkdir -p "${target_dir}/counted"
cd "${target_dir}" || exit

## Count
for sample_name in "${sample_names[@]}"; do
    cellranger count --id="${sample_name}" --sample="${sample_name}" --fastqs="$target_dir/raw" \
        --transcriptome="$reference_genome" --no-bam
done

cd "${target_dir}/raw" || exit
## Aggregate
for sample_name in "${sample_names[@]}"; do
    cellranger aggr --id="${sample_name}" --csv=meta/aggregation_"${sample_name}".csv
done
