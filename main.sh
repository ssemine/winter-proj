#!/bin/bash
# main.sh
# -------
# Main script that runs transform.sh and run.sh for each gene in the input file.
# ------------------------------------------------------------------------------
# Usage: ./main.sh --infile <input_file> --bfile <bed_file> --maf <maf> --chr <chr_num> \
# 			[ --pval <p_value> | --genes <gene_list> | --snps <snp_list> | --log <log_file> | --gene_dir <gene_dir> ]

source definitions/constants.sh
source definitions/file_indices.sh
source definitions/log_messages.sh

gene_list="$GENE_LIST"
gene_dir="$GENE_DIR"
snp_dir="$SNP_DIR"
genes="$GENES_ALL" 
snps="$SNPS_ALL"
p_val="$P_VALUE_THRESHOLD"

log_dir="$LOG_DIR"
log_file="$LOG_FILE"
summary_file="$SUMMARY_FILE"
snp_csv_header="$SNP_CSV_HEADER"
snp_count_file="$SNP_COUNT_FILE"

source definitions/functions.sh
mkdir -p "$log_dir"

 
while [[ $# -gt 0 ]]; do
    case $1 in
        --infile)
            infile="$2"
            shift 2
            ;;
        --bfile)
            bfile="$2"
            shift 2
            ;;
        --maf)
            maf="$2"
            shift 2
            ;;
        --pval)
            p_val="$2"
            shift 2
            ;;
		--chr)
			chr="$2"
			shift 2
			;;
        --genes)
      	    genes="$2"
            shift 2
            ;;
		--snps)
	    	snps="$2"
	    	shift 2
	    	;;
		--log)
			log_file="$2"
			shift 2
			;;
		--gene_dir)
			gene_dir="$2"
			shift 2
			;;
        *)
            echo "$ERROR_INVALID_ARGUMENT $1"
            exit 1
            ;;
    esac
done


summary_file="${log_file%.log}_summary.log"
summary_log() {
	local message="$1"
	echo "$message" >> "$log_dir/$summary_file"
}



touch "$log_dir/$log_file"
touch "$log_dir/$summary_file"

log "$LOG_WELCOME_MESSAGE"
log_lines 2

if [[ -z "$chr" ]]; then
    { log "$ERROR_CHR_NUM_NOT_PROVIDED"; exit 1; }
fi

log "$(printf "$LOG_PARAMETERS" "$infile" "$bfile" "$maf" "$p_val" "$chr" "$genes" "$snps" "$log_file")"
log_lines 1


if [ "$genes" == "$GENES_ALL" ] && [ ! -f "$genes" ]; then
	awk -v gidx="$INPUT_GENE_NAME_IDX" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list" \
		|| { log "$ERROR_GENE_LIST"; exit 1; }
else
	cat "$genes" > "$gene_list"
fi
summary_log "Number of genes: $(wc -l < $gene_list)"
echo "" >> "$log_dir/$summary_file"

touch "$snp_count_file"
echo "$snp_csv_header" > "$snp_count_file"
log "$LOG_CREATED_FILES $snp_count_file"

if [ -f "$snps" ]; then
    if ! awk -v sidx="$INPUT_SNP_ID_IDX" 'NR==FNR {snps[$1]; next} $sidx in snps' "$snps" "$infile" | sort | uniq -c | awk '{ print $2, $1 }' >> "$snp_count_file"; then
        { log "$ERROR_FILTERED_SNP_COUNT_FILE"; exit 1; }
    fi
else
    if ! awk -v sidx="$INPUT_SNP_ID_IDX" '{ print $sidx }' "$infile" | sort | uniq -c | awk '{ print $2, $1 }' >> "$snp_count_file"; then
        { log "$ERROR_SNP_COUNT_FILE"; exit 1; }
    fi
fi

mkdir -p "$gene_dir"
mkdir -p "$snp_dir"

while IFS= read -r line; do
	log_lines 1
	log "$LOG_CALLING_TRANSFORM $line"
	./transform.sh "$line" \
		"$infile" \
		"$gene_dir" \
		"$snps" \
		"$chr" \
		"$log_file" \
		"$log_dir" \
		"$bfile" \
		"input"
	log "$LOG_MA_TRANSFORMED $line"
	log_lines 1
	log "$LOG_CALLING_RUN $line"
	summary_log "Gene: $line"
	./run.sh "$line" \
		"$bfile" \
		"$chr" \
		"$maf" \
		1 \
		"$p_val" \
		"$log_file" \
		"$log_dir" \
		"$gene_dir" \
		"$snp_dir" \
		"$summary_file" \
		"$snps"
	log "$LOG_RUN_FINISHED $line"
	echo "" >> "$log_dir/$summary_file"
done < "$gene_list" || { log "$ERROR_READ_GENE_LIST"; exit 1; }
log_lines 1
log "$LOG_END_MESSAGE"