#!/bin/bash

# run.sh
# ------
# Fetches lowest p-value SNP from .ma file and runs GCTA conditional analysis. 
# Calls transform.sh to transform the GCTA's ouput .cma.cojo to .ma file, passes it to run.sh until all SNPs are with p-value < threshold are fetched.
# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Usage: ./run.sh gene_name bfile chr maf idx p_val log_file log_dir gene_dir snp_dir summary_file
# ------------------------------------------------------------------------------------------------

source definitions/constants.sh
source definitions/file_indices.sh
source definitions/log_messages.sh

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
outfile=$(printf "$GCTA_OUTFILE_NAME" "$gene_name" "$idx")
infile=$(printf "$MA_FILE_NAME" "$gene_name")

source definitions/functions.sh

log "$LOG_STARTING_RUN $gene_name"
log "$LOG_ITERATION_NUM $idx"

# If idx = 1, it means .ma file is used to fetch the lowest p-value
if [ $idx -eq 1 ]; then
	read_file="$infile"
else
	read_file=$(printf "$MA_FILE_NAME_IDX" "$gene_name" $prev_idx)
fi
log "$LOG_READING_FROM $gene_dir/$read_file"

top_snp_file=$(printf "$TOP_SNP_FILE" "$snp_dir" "$gene_name" "$idx")
touch "$top_snp_file"

awk -v col="$MA_P_VALUE_IDX" \
    -v id_col="$MA_SNP_ID_IDX" \
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
    || { log "$ERROR_AWK_WRITE $top_snp_file"; exit 1; }

# Checks if top snp file is empty
has_snp=$(wc -l < "$top_snp_file")

if [ "$has_snp" -eq 1 ]; then
	log "$(printf $LOG_TOP_SNP $gene_name $(cat $top_snp_file))"
	./gcta64 --bfile "$bfile" \
        --chr "$chr" \
        --maf "$maf" \
        --cojo-file "$gene_dir/$read_file" \
		--cojo-cond "$top_snp_file" \
        --out "$gene_dir/$outfile"
    ./transform.sh "$gene_name" \
        "$(printf $TRANSFORM_CMA_FILE_NAME $gene_dir $outfile)" \
        "$gene_dir" \
        "$snps" \
        "$chr" \
        "$log_file" \
        "$log_dir" \
        "$bfile" \
        "$CMA_IDENTIFIER" \
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
	log $(printf "$LOG_TOTAL_SNPS" "$gene_name" "$prev_idx")
    summary_log $(printf "$LOG_TOTAL_SNPS" "$gene_name" "$prev_idx")
fi