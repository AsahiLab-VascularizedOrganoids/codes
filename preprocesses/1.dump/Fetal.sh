#!/bin/bash

####### Params ########################
dump_target_dir="$HOME/resources/Fetal/raw/"
#######################################

mkdir -p "${dump_target_dir}"
cd "${dump_target_dir}" || exit
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162170/suppl/GSE162170_rna_counts.tsv.gz" -O "count.tsv"
