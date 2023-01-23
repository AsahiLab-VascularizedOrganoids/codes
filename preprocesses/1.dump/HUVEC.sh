#!/bin/bash

####### Params ########################
dump_target_dir="$HOME/resources/HUVEC/raw/"
srrIDs=(
    "SRR9047768" "SRR10970638" "SRR10970639" #HUVEC_hCO_d65
    "SRR9047769" "SRR10970640" "SRR10970641" #HUVEC_vhCO_d65
    "SRR9047770" "SRR10970642" "SRR10970643" #HUVEC_hCO_d100
    "SRR9047771" "SRR10970644" "SRR10970645" #HUVEC_vhCO_d100
)
#######################################

mkdir -p "${dump_target_dir}"
cd "${dump_target_dir}" || exit

for srrID in "${srrIDs[@]}"; do
    echo "Starting dump ${srrID}"
    fastq-dump --split-files --origfmt --gzip -O "${dump_target_dir}" "${srrID}"
done
