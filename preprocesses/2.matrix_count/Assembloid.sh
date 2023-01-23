#!/bin/bash

####### Params ########################
sample_name="VOr"
target_dir="$HOME/resources/$sample_name/"
reference_genome="${HOME}/ref/refdata-gex-GRCh38-2020-A"
#######################################

cd "${target_dir}" || exit

## Count
bamtofastq --nthreads=4 "raw" "fq"
## Aggregate
cellranger count --id="$sample_name" --sample="$sample_name" --fastqs="fq" \
    --transcriptome="$reference_genome" --no-bam
