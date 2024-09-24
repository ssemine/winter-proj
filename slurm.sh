#!/bin/bash

#SBATCH --account=CHANGE_ACCOUNT_GROUP_NAME
#SBATCH --mem=1000GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=10
#SBATCH --partition=general
#SBATCH --time=06:00:00
#SBATCH --job-name=scojo

infile="$1"
bfile="$2"
maf="$3"
p_val="$4"
chr="$5"
exclude_qtl_type="$6"
run_dir="$7"

./main.sh --infile "$infile" \
    --bfile "$bfile" \
    --maf "$maf" \
    --pval "$p_val" \
    --chr "$chr" \
    --exclude_qtl_type "$exclude_qtl_type" \
    --run_dir "$run_dir" \
    --log log_file \
    --genes all
