#!/bin/bash

####### Params ########################
dump_target_dir="$HOME/resources/VOr/raw/"
#######################################

mkdir -p "${dump_target_dir}"
cd "${dump_target_dir}" || exit
wget "https://sra-pub-src-1.s3.amazonaws.com/SRR15992286/V.bam.1" -O "VOR1.bam"
wget "https://sra-pub-src-2.s3.amazonaws.com/SRR15992285/V2.bam.1" -O "VOR2.bam"
