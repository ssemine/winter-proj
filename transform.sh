#!/bin/bash

# DEBUG

# This script transforms input file to .ma file. 

# CHECK INFILE -> IT SHOULD BE IN GENE_DIR - fixed
# NEED TO ENSURE .ma is transformed


# Command line arguments indices
gene_name="$1"
infile="$2"
gene_dir="$3"
snp_id_idx="$4"
chr_idx="$5"
allele_one_idx="$6"
allele_two_idx="$7"
freq_idx="$8"
gene_name_idx="$9"
effect_size_idx="${10}"
se_idx="${11}"
p_value_idx="${12}"
snps="${13}"
chr_num="${14}"
log_file="${15}"
log_dir="${16}"
bfile="${17}"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_dir/$log_file"
}

log "Starting transform.sh for $gene_name" 


# Variables declarations
name="$gene_name.ma"
name_final="${gene_name}_input.ma"
columns="SNP A1 A2 freq b se p N"

# Creates temporary and final .ma files

if ! touch "$gene_dir/$name"; then
    log "transform.sh Error: touch could not create $gene_dir/$name"
    exit 1
fi

if ! touch "$gene_dir/$name_final"; then
    log "transform.sh Error: touch could not create $gene_dir/$name_final"
    exit 1
fi

log "    .ma indices:"
log "        SNP ID: $snp_id_idx"
log "        Chromosome: $chr_idx"
log "        Allele 1: $allele_one_idx"
log "        Allele 2: $allele_two_idx"
log "        Frequency: $freq_idx"
log "        Gene Name: $gene_name_idx"
log "        Effect Size: $effect_size_idx"
log "        Standard Error: $se_idx"
log "        P-Value: $p_value_idx"

# Writes required information from input column to .ma file if gene name matches
log "	Writing data to $gene_dir/$name"
if [ -f "$snps" ]; then
    log "	Using $snps to filter SNPs"
    awk -F' ' -v col1="$snp_id_idx" \
        -v col2="$allele_one_idx" \
        -v col3="$allele_two_idx" \
        -v col4="$freq_idx" \
        -v col5="$effect_size_idx" \
        -v col6="$se_idx" \
        -v col7="$p_value_idx" \
        -v gidx="$gene_name_idx" \
        -v cidx="$chr_idx" \
        -v gene="$gene_name" \
        -v gene_dir="$gene_dir" \
        -v name="$name" \
        -v snps="$snps" \
        -v chr_num="$chr_num" \
        'BEGIN {
            while ((getline < snps) > 0) {
                snp_list[$1]
            }
        }
        {
            if ($gidx == gene && ($col1 in snp_list) && $cidx == chr_num) {
                print $col1, $col2, $col3, $col4, $col5, $col6, $col7 >> (gene_dir "/" name)
            }
        }' "$infile" \
        || { log "transform.sh Error: awk could not write data to $gene_dir/$name"; exit 1; }
else
    log "	No SNP filter"
    awk -F' ' -v col1="$snp_id_idx" \
        -v col2="$allele_one_idx" \
        -v col3="$allele_two_idx" \
        -v col4="$freq_idx" \
        -v col5="$effect_size_idx" \
        -v col6="$se_idx" \
        -v col7="$p_value_idx" \
        -v gidx="$gene_name_idx" \
        -v cidx="$chr_idx" \
        -v gene="$gene_name" \
        -v gene_dir="$gene_dir" \
        -v name="$name" \
        -v name_chr="${gene_name}_chr.txt" \
        '{
            if ($gidx == gene) {
                print $col1, $col2, $col3, $col4, $col5, $col6, $col7 >> (gene_dir "/" name)
                print gene_name, $cidx >> (gene_dir "/" name_chr)
            }
        }' "$infile" \
        || { log "transform.sh Error: awk could not write data to $gene_dir/$name"; exit 1; }
fi

log "	Data written to $gene_dir/$name"
#log "   Inspecting data before sorting"
#head -n 5 "$gene_dir/$name" >> "$log_dir/$log_file"

# Sorts the temporary file by SNP
sort -k 1 "$gene_dir/$name" -o "$gene_dir/$name" \
    || { log "transform.sh Error: sort could not sort $gene_dir/$name"; exit 1; }

sample_size=$(wc -l < "$bfile.fam")

# Adds the sample_size number to each row, writes data to a new file
awk -v sample_size="$sample_size" '{print $0, sample_size}' "$gene_dir/$name" > "$gene_dir/$name_final" \
    || { log "transform.sh Error: awk could not write data to $gene_dir/$name_final"; exit 1; }

log "	Sample size $sample_size added to $gene_dir/$name_final"


# Removes temporary .ma file
# rm "$gene_dir/$name"

# Adds column headers to .ma file
awk -v "cols=$columns" 'BEGIN{print cols}1' "$gene_dir/$name_final" > temp \
	|| { log "transform.sh Error: awk could not add columns to $gene_dir/$name_final"; exit 1; }
	
mv temp "$gene_dir/$name_final"
log "	Columns added to $gene_dir/$name_final"