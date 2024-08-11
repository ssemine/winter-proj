#!/bin/bash

# POLISH
# Usage:
# 	main.sh --infile infile --bfile bfile --maf maf --chr chromosome_number [ --p-value p_val | --genes gene_names | --snps snp_names | --log log_file | --gene_dir gene_directory ]

# Arguments: 
# 	infile: input file
# 	bfile: bed files
# 	maf: minor allele frequency
# 	chr: chromosome number
#	p_val: p-value threshold (optional)
# 	genes: gene names (optional)
# 	snps: snp names (optional)
# 	log: log file (optional)
#	gene_dir: gene directory (optional)

source definitions/file_indices.sh

# CONSTANTS
snp_id_idx=1
chr_idx=2
allele_one_idx=4
allele_two_idx=5
freq_idx=6
gene_name_idx=7
effect_size_idx=12
se_idx=13
p_value_idx=14
gene_list="gene_list.txt"
gene_dir="cojo_files"
snp_dir="snp_files"
genes="all" 
snps="all"
log_dir="logs"
p_val="5e-12"
snp_csv_header="SNP,Count"
snp_count_file="snp_count.csv"
summary_file="summary.log"


mkdir -p "$log_dir"

log_file="$(date '+%Y-%m-%d %H:%M:%S').log"



# Function to log messages
# Usage: log "message"
log() {
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo "$file_name:$line_number - $message" >> "$log_dir/$log_file"
}
log_lines() {
    local num_lines="$1"
    for ((i = 0; i < num_lines; i++)); do
        echo "" >> "$log_dir/$log_file"
    done
}

# Allows for dynamic argument assignment
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
            echo "Error: invalid argument: $1"
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

log "Start of script"
log_lines 2

if [[ -z "$chr" ]]; then
    { log "Error: chromosome number not provided"; exit 1; }
fi

log "Parameters selected: "
log "Input file: $infile"
log "BED files: $bfile"
log "Minor allele frequency: $maf"
log "p-value threshold: $p_val"
log "Genes: $genes"
log "SNPs: $snps"
log "Log: $log_file"
log "Opened $infile"
log_lines 1
# Gene selection
if [ "$genes" == "all" ] && [ ! -f "$genes" ]; then
	awk -v gidx="$INPUT_GENE_NAME_IDX" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list" \
		|| { log "Error: unable to create gene list"; exit 1; }
else
	cat "$genes" > "$gene_list"
fi
summary_log "Number of genes: $(wc -l < $gene_list)"
echo "" >> "$log_dir/$summary_file"



# SNP files
touch "$snp_count_file"
echo "$snp_csv_header" > "$snp_count_file"
# Check if $snps is a file and filter SNPs accordingly
if [ -f "$snps" ]; then
    if ! awk -v sidx="$INPUT_SNP_ID_IDX" 'NR==FNR {snps[$1]; next} $sidx in snps' "$snps" "$infile" | sort | uniq -c | awk '{ print $2, $1 }' >> "$snp_count_file"; then
        { log "Error: unable to create filtered snp count file"; exit 1; }
    fi
else
    if ! awk -v sidx="$INPUT_SNP_ID_IDX" '{ print $sidx }' "$infile" | sort | uniq -c | awk '{ print $2, $1 }' >> "$snp_count_file"; then
        { log "Error: unable to create snp count file"; exit 1; }
    fi
fi

log "Created file snp_count.txt"

mkdir -p "$gene_dir"
mkdir -p "$snp_dir"

while IFS= read -r line; do
	log_lines 1
	log "Calling transform.sh for $line"
	./transform.sh "$line" \
		"$infile" \
		"$gene_dir" \
		"$snps" \
		"$chr" \
		"$log_file" \
		"$log_dir" \
		"$bfile" \
		"input"
	log ".ma for $line transformed"
	log_lines 1
	log "Calling run.sh for $line"
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
	log "run.sh for $line finished"
	echo "" >> "$log_dir/$summary_file"
done < "$gene_list" || { log "Error: unable to read gene list"; exit 1; }
log_lines 1
log "End of script"