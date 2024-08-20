#!/bin/bash

# test
# main.sh
# -------
# Main script that runs transform.sh and run.sh for each gene in the input file.
# ------------------------------------------------------------------------------
# Usage: ./main.sh --infile <input_file> --bfile <bed_file> --maf <maf> --chr <chr_num> \
# 			[ --pval <p_value> | --genes <gene_list> | --snps <snp_list> | --log <log_file> | --gene_dir <gene_dir> ]
# -------------------------------------------------------------------------------------------------------------------

source definitions/constants.sh
source definitions/file_indices.sh
source definitions/log_messages.sh

gene_list="$GENE_LIST"
gene_dir="$GENE_DIR"
snp_dir="$SNP_DIR"
genes="$GENES_ALL" 
snps="$SNPS_ALL"
p_val="$P_VALUE_THRESHOLD_PER_CHR"

log_dir="$LOG_DIR"
log_file="$LOG_FILE"
summary_file="$SUMMARY_FILE"
snp_csv_header="$SNP_CSV_HEADER"
snp_count_file="$SNP_COUNT_FILE"

run_dir="$RUN_DIR"

results_file="$RESULTS_FILE_NAME"


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
		--run_dir)
			run_dir="$2"
			shift 2
			;;
		--clean)
			clean="$2"
			shift 2
			;;
        *)
            echo "$ERROR_INVALID_ARGUMENT $1"
            exit 1
            ;;
    esac
done

if [[ -d "$run_dir" ]]; then
    echo "$(printf "$RUN_DIR_EXISTS" "$run_dir")"
	exit 1
fi

mkdir -p "$run_dir"
summary_file="${log_file%.log}_summary.log"
source definitions/functions.sh
cd "$run_dir"
mkdir -p "$log_dir"
touch "$log_dir/$log_file"
touch "$log_dir/$summary_file"
touch "$results_file"
echo "$RESULTS_FILE_HEADER" > "$results_file"
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

mkdir -p "$gene_dir"
mkdir -p "$snp_dir"

if [ "$p_val" = "$P_VALUE_THRESHOLD_PER_CHR" ]; then
	num_snps=$(awk -v snp_col="$INPUT_SNP_ID_IDX" '{ print $snp_col }' "$infile" | sort | uniq | wc -l)
	p_val=$(echo "scale=$P_VALUE_PRECISION; $P_VALUE_NUMERATOR / $num_snps" | bc)
	p_val=$(printf "%.${P_VALUE_PRECISION}f\n" "$p_val")
fi

while IFS= read -r line; do
	log_lines 1
	log "$LOG_CALLING_TRANSFORM $line"
	"$PATH_TO_TRANSFORM_SH" "$line" \
		"$infile" \
		"$gene_dir" \
		"$chr" \
		"$log_file" \
		"$log_dir" \
		"$bfile" \
		"$INPUT_IDENTIFIER" \
		"$PATH_TO_DEFINITIONS" \
		"$snps" \
		|| { log "$ERROR_TRANSFORM $line"; exit 1; }
	ma_file_reference="$(printf "$MA_FILE_NAME_REFERENCE" "$gene_dir/$line")"
	cp "$gene_dir/$line.ma" "$ma_file_reference"
	log "$LOG_MA_TRANSFORMED $line"

	# Calculate p-value per gene (consider to remove)
	if [ p_val = "$P_VALUE_THRESHOLD_PER_GENE" ]; then
		num_snps=$(wc -l < "$gene_dir/$line.ma")
		p_val=$(echo "scale=$P_VALUE_PRECISION; $P_VALUE_NUMERATOR / $num_snps" | bc)
		p_val=$(printf "%.${P_VALUE_PRECISION}f\n" "$p_val")
	fi
	summary_log "p-value for $line $p_val"
	log_lines 1
	log "$LOG_CALLING_RUN $line"
	summary_log "Gene: $line"
	"$PATH_TO_RUN_SH" "$line" \
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
		"$PATH_TO_DEFINITIONS" \
		"$results_file" \
		|| { log "$ERROR_RUN_FAILED $line"; exit 1; }

	# Join GCTA logs by gene and run.sh idx
	gcta_log_file="$(printf "$GCTA_LOG_FILE" "$line")"
	touch "$log_dir/$gcta_log_file"
	cat "$gene_dir/$line"*.log > "$log_dir/$gcta_log_file"
	rm "$gene_dir/$line"*.log

	log "$LOG_RUN_FINISHED $line"

	# Clean up
	rm "$ma_file_reference"
	rm *".$SNP_LINE_EXTENTION"

	echo "" >> "$log_dir/$summary_file"
done < "$gene_list" || { log "$ERROR_READ_GENE_LIST"; exit 1; }


# Clean up
if [ "$clean" = "all" ]; then 
	rm "$gene_list"
	rm -rf "$gene_dir"
	rm -rf "$snp_dir"
fi

log_lines 1
log "$LOG_END_MESSAGE"