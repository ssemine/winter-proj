#!/bin/bash

#SBATCH --account=INSRT
#SBATCH --mem=1000GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=10
#SBATCH --partition=general
#SBATCH --time=03:00:00
#SBATCH --job-name=testcojo

infile="$1"
bfile="$2"
maf="$3"
p_val="$4"
chr="$5"
exclude_qtl_type="$6"

./main.sh --infile "$infile" \
    --bfile "$bfile" \
    --maf "$maf" \
    --pval "$p_val" \
    --chr "$chr" \
    --exclude_qtl_type "$exclude_qtl_type" \
    --run_dir "$bfile" \
    --log log_file \
    --genes all