#!/bin/bash

infile="$1"
bfiles_dir="$2"
maf="$3"
p_val="$4"
exclude_qtl_type="$5"

bfiles=$(ls "$chr_file_dir" | cut -d. -f1 | sort | uniq)


for bfile in "$bfiles"; do
    chr=$(echo "$bfile" | grep -oP '(?<=chrom)\d+')
    sbatch slurm.sh "$infile" \
        "$bfile" \
        "$maf" \
        "$p_val" \
        "$chr" \
        "$exclude_qtl_type"