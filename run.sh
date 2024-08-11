#!/bin/bash

# POLISH

# $1 gene $2 - bfile, $3 - chr, $4 - maf, $5 - run index, $6 p-value threshold
# outputs .cma file
#

gene_name="$1"
bfile="$2"
chr="$3"
maf="$4"
idx="$5"
p_val="$6"
log_file="$7"
log_dir="$8"
gene_dir="$9"
snp_dir="${10}"
summary_file="${11}"
prev_idx="$((idx - 1))"
next_idx="$((idx + 1))"
outfile=$(printf "%s_%s" "$gene_name" "$idx")
infile=$(printf "%s.ma" "$gene_name")


log() {
    local gene="$gene_name"
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo "$file_name:$line_number - $gene: $message" >> "$log_dir/$log_file"
}
log_lines() {
    local num_lines="$1"
    for ((i = 0; i < num_lines; i++)); do
        echo "" >> "$log_dir/$log_file"
    done
}
summary_log() {
	local message="$1"
	echo "$message" >> "$log_dir/$summary_file"

}
log "Iteration number: $idx"

# If idx = 1, it means .ma file is used to fetch the lowest p-value
if [ $idx -eq 1 ]; then
	read_file="$infile"
else
	read_file=$(printf "%s_%s.ma" "$gene_name" $prev_idx)
fi
snp_col=1
p_col=7
log "Reading from $gene_dir/$read_file"

top_snp_file=$(printf "%s/%s_%s.snplist" "$snp_dir" "$gene_name" "$idx")
touch "$top_snp_file"

awk -v col="$p_col" \
    -v id_col="$snp_col" \
    -v thresh="$p_val" \
    -v snp_file="$top_snp_file" \
    -v summary_file="$log_dir/$summary_file" \
    'NR > 1 && $col < thresh { 
        if (min == "" || $col < min) { 
            min = $col; 
            id = $id_col 
        } 
    } 
    END { 
        if (min != "" && min < thresh) 
            print id > snp_file
            print id, min >> summary_file
    }' "$gene_dir/$read_file" \
    || { log "Error: awk unable to create $top_snp_file"; exit 1; }

# Checks if top snp file is empty
has_snp=$(wc -l < "$top_snp_file")

if [ "$has_snp" -eq 1 ]; then
	log "Top SNP for $gene_name: $(cat $top_snp_file)"
	./gcta64 --bfile "$bfile" \
        --chr "$chr" \
        --maf "$maf" \
        --cojo-file "$gene_dir/$read_file" \
		--cojo-cond "$top_snp_file" \
        --out "$gene_dir/$outfile"
    ./transform.sh "$gene_name" \
        "$gene_dir/$outfile.cma.cojo" \
        "$gene_dir" \
        "$snps" \
        "$chr" \
        "$log_file" \
        "$log_dir" \
        "$bfile" \
        "cma" \
        "$gene_dir/$infile" \
        "$idx"
	./run.sh "$gene_name" \
		"$bfile" \
		"$chr" \
		"$maf" \
		"$next_idx" \
		"$p_val" \
		"$log_file" \
		"$log_dir" \
		"$gene_dir" \
        "$snp_dir" \
        "$summary_file"
else
	log "Total SNPs for $gene_name: $prev_idx"
    summary_log "Total SNPs for $gene_name: $prev_idx"
fi