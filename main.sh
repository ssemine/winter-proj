#!/bin/bash

# main.sh
# -------
# Main script that runs transform.sh and run.sh for each gene in the input file.
# ------------------------------------------------------------------------------
# Usage: 				./main.sh \
# Required arguments:		--infile <input_file> \
# 							--bfile <bed_file> \
#							--maf <maf> \
#							--chr <chr_num> \
# Optional arguments:		--pval <p_value_method> or <p_value> \
#							--genes <gene_list> \
#							--snps <snp_list> \
#							--log <log_file> \
#							--gene_dir <gene_dir> \
#							--run_dir <run_dir> \
#							--exclude_qtl_type <qtl_type>
# ------------------------------------------------------------------------------
# NOTE: Current implementation assumes running with no trans eQTLs.

# Variables are sourced
source definitions/constants.sh
source definitions/file_indices.sh
source definitions/log_messages.sh

# Local variables are set
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

# Command line arguments are parsed
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
		--exclude_qtl_type)
			exclude_qtl_type="$2"
			shift 2
			;;
        *)
            echo "$ERROR_INVALID_ARGUMENT $1"
            exit 1
            ;;
    esac
done

# Checks if the run directory already exists, if so, exits
if [[ -d "$run_dir" ]]; then
    echo "$(printf "$RUN_DIR_EXISTS" "$run_dir")"
	exit 1
fi

# Checks if chromosome number is provided
if [[ -z "$chr" ]]; then
    { log "$ERROR_CHR_NUM_NOT_PROVIDED"; exit 1; }
fi

# Creates necessary directories and files
mkdir -p "$run_dir"
summary_file="${log_file%.log}_summary.log"
source definitions/functions.sh
cd "$run_dir"
mkdir -p "$log_dir"
mkdir -p "$log_dir/$GCTA_LOG_DIR"
touch "$log_dir/$log_file"
touch "$log_dir/$summary_file"
touch "$results_file"
echo "$RESULTS_FILE_HEADER" > "$results_file"

# Logs the start of the script
log "$LOG_WELCOME_MESSAGE"
log_lines 2
log "$(printf "$LOG_PARAMETERS" "$infile" "$bfile" "$maf" "$p_val" "$chr" "$genes" "$snps" "$log_file")"
log_lines 1

# Creates gene_list file
if [ "$genes" == "$GENES_ALL" ] && [ ! -f "$genes" ]; then
	awk -v gidx="$INPUT_GENE_NAME_IDX" \
		-v cidx="$INPUT_CHR_IDX" \
		-v chr="$chr" \
		'{
			if ($cidx == chr) {
				print $gidx
			}
		}' "$infile" | sort -k 1 | uniq > "$gene_list" \
		|| { log "$ERROR_GENE_LIST"; exit 1; }
else
	cat "$genes" > "$gene_list"
fi

# Excludes pares QTL type
awk -v exclude_qtl_type="$exclude_qtl_type" \
	-v chr="$chr" \
	-v chr_idx="$INPUT_CHR_IDX" \
	-v qtl_type_idx="$INPUT_QTL_TYPE_IDX" \
	'{
		if (($chr_idx == chr) && ($qtl_type_idx != exclude_qtl_type)) {
			print $0
		}
	}' "$infile" > "$INFILE_COPY"
infile="$INFILE_COPY"

# If real chromosome number is outside the GCTA range, change it to 1
original_chr="$chr"
if [[ "$chr" -ge 23]]; then
	chr="$SUBSTITUTE_CHR"
	awk -v chr_idx="$INPUT_CHR_IDX" \
    		-v snp_idx="$INPUT_SNP_ID_IDX" \
    		-v chr="$chr" \
    		'{
        		split($snp_idx, snp, ":");
        		$(chr_idx) = chr;
        		$(snp_idx) = chr ":" snp[2];
        		print $0;
    		}' "$infile" > "$infile.tmp"
	rm "$infile"
	mv "$infile.tmp" "$infile"
fi

# Create SNP helper list 
awk -v snp="$INPUT_SNP_ID_IDX" \
	-v pos_snp="$INPUT_POS_SNP_IDX" \
	-v pos_gene="$INPUT_POS_GENE_IDX" \
	-v gene_name="$INPUT_GENE_NAME_IDX" \
	-v strand="$INPUT_STRAND_IDX" \
	-v qtl_type="$INPUT_QTL_TYPE_IDX" \
	'{
		print $snp, $pos_snp, $pos_gene, $gene_name, $strand, $qtl_type
	}' "$infile" | sort -k 1 | uniq > "$SNP_HELPER_LIST"

summary_log "Number of genes: $(wc -l < $gene_list)"
echo "" >> "$log_dir/$summary_file"

# Create gene and snp directories
mkdir -p "$gene_dir"
mkdir -p "$snp_dir"

# Compute p-value threshild per chromosome
if [ "$p_val" = "$P_VALUE_THRESHOLD_PER_CHR" ]; then
	num_snps=$(awk -v snp="$INPUT_SNP_ID_IDX" \
		-v cidx="$INPUT_CHR_IDX" \
		-v chr="$chr" \
		'{
			if ($cidx == chr) {
				print $snp
			}
		}' "$infile" | sort -k 1 | uniq | wc -l)
	p_val=$(echo "scale=$P_VALUE_PRECISION; $P_VALUE_NUMERATOR / $num_snps" | bc)
	p_val=$(printf "%.${P_VALUE_PRECISION}f\n" "$p_val")
fi

# Main loop
while IFS= read -r line; do
	log_lines 1
	log "$LOG_CALLING_TRANSFORM $line"

	# Transform the input data into .ma file format
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

	# Copy .ma file to reference file, used in run.sh
	ma_file_reference="$(printf "$MA_FILE_NAME_REFERENCE" "$gene_dir/$line")"
	cp "$gene_dir/$line.ma" "$ma_file_reference"
	log "$LOG_MA_TRANSFORMED $line"

	# Calculate p-value per gene (consider to remove)
	if [ p_val = "$P_VALUE_THRESHOLD_PER_GENE" ]; then
		num_snps=$(cat "$gene_dir/$line.ma" | awk -v col="$MA_SNP_ID_IDX" '{ print $col }' | sort -k 1 | uniq | wc -l)
		p_val=$(echo "scale=$P_VALUE_PRECISION; $P_VALUE_NUMERATOR / $num_snps" | bc)
		p_val=$(printf "%.${P_VALUE_PRECISION}f\n" "$p_val")
	fi

	summary_log "p-value for $line $p_val"
	log_lines 1
	log "$LOG_CALLING_RUN $line"
	summary_log "Gene: $line"

	# Run GCTA conditional analysis for the gene
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
	gcta_log_file="$log_dir/$GCTA_LOG_DIR/$gcta_log_file"
	touch "$gcta_log_file"
	cat "$gene_dir/$line"*.log > "$gcta_log_file"
	rm "$gene_dir/$line"*.log

	log "$LOG_RUN_FINISHED $line"

	# Clean up
	rm "$ma_file_reference"
	rm *".$SNP_LINE_EXTENTION"

	echo "" >> "$log_dir/$summary_file"
done < "$gene_list" || { log "$ERROR_READ_GENE_LIST"; exit 1; }

# Fix chromosome and SNP IDs in results + skipping header
chr="$original_chr"
awk -v chr="$chr" \
    -v chr_idx="$OUTPUT_CHR_IDX" \
    -v snp_idx="$OUTPUT_SNP_ID_IDX" \
    'NR == 1 { print $0; next } 
    {
        split($snp_idx, snp, ":");
        $(chr_idx) = chr;
        $(snp_idx) = chr ":" snp[2];
        print $0;
    }' "$results_file" > "$results_file.tmp"
rm "$results_file"
mv "$results_file.tmp" "$results_file"

# Remove tempotaory bfiles
if [[ "$chr" -ge "$GCTA_MAX_CHR" ]]; then
	rm "$bfile.bim"
	mv "$bfile.bim.tmp" "$bfile.bim"
fi

log_lines 1
log "$LOG_END_MESSAGE"
