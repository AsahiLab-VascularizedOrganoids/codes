#!/bin/bash

####### Params ########################
dump_target_dir="$HOME/resources/hETV2/raw/"
#######################################

mkdir -p "${dump_target_dir}"
cd "${dump_target_dir}" || exit
wget "https://sra-pub-src-2.s3.amazonaws.com/SRR9661365/Sample1.bam.1" -O "hETV2_hCO.bam"
wget "https://sra-pub-src-2.s3.amazonaws.com/SRR9661366/Sample2.bam.1" -O "hETV2_vhCO.bam"
