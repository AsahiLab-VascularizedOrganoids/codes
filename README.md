# Vascularization_organoid > Meta Analysis

## About This

Project directory of
sato et al., Comparative Single-Cell Analysis of Vascularization for Human Cortical Organoids

## Analysis Flow

1. resources/1.fq
   dump fastq file(s) from SRA database

2. [optional]resources/2.cellranger
   run cellranger to generate single-cell matrix

3. resources/3.mat
   generate matrix from cellranger output > pre-process.py

4. resources/4.integration
   integrate all matrix using `Seurat`

5. integrated
   analyze single cell data using `Seurat`, `SingleCellPipeline`, `dplyr`, `slingshot` etc...

## Samples

hETV2, HUVEC, assembloid (VOr) and fetals.
The following fetal samples are included on "fetals".

- fetalPCW16_GSE162170
- fetalPCW20_GSE162170
- fetalPCW21_GSE162170
- fetalPCW24_GSE162170

## Methods

1. SCT: True

2. Feature_Cutoff: 1000

3. Resolution: 1.2, 2.0

## Directories & Files

## Requirements
