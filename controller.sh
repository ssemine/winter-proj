#!/bin/bash

# controller.sh
# -------------

infile=$1
bfile=$2
maf=$3
chr=$4
genes=$5
log_file=$6
p_val=$7
run_dir=$8

source definitions/constants.sh


split -l $GENES_PER_MAIN -d --additional-suffix=.genelist "$genes" "${genes}_part_"

# Loop through the generated files and submit a job for each one
for file in ${genes}_part_*.genelist; 
do
    sbatch slurm.sh "$infile" "$bfile" "$maf" "$chr" "$file" "$log_file" "$p_val" "$run_dir.$file"
done



