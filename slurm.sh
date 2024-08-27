#!/bin/bash
#SBATCH --account=insert_accounting_group_here
#SBATCH --mem=1000GB
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=10
#SBATCH --time=01:00:00
#SBATCH --job-name=gctacojo

infile=$1
bfile=$2
maf=$3
chr=$4
genes=$5
log_file=$6
p_val=$7
run_dir=$8

./main.sh --infile "$infile" \
        --bfile "$bfile" \
        --maf "$maf" \
        --chr "$chr" \
        --genes "$genes" \
        --log "$log_file" \
        --pval "$p_val" \
        --run_dir "$run_dir"