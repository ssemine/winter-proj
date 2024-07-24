#!/bin/bash

# DEBUG

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

# CONSTANTS
snp_id_idx=1
chr_idx=2
allele_one_idx=4
allele_two_idx=5
freq_idx=6
gene_name_idx=7
p_value_idx=14
gene_list="gene_list.txt"
gene_dir="cojo_files"
genes="all" 
snps="all"
log_dir="logs"
p_val="5e-12"

mkdir -p "$log_dir"

log_file="$log_dir/$(date '+%Y-%m-%d %H:%M:%S').log"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_dir/$log_file"
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

touch "$log_dir/$log_file"

if [[ -z "$chr" ]]; then
    log "main.sh Error: chromosome number not provided" && exit 1
fi

log "Running main.sh..."
log "	Parameters selected: "
log "		Input file: $infile"
log "		BED files: $bfile"
log "		Minor allele frequency: $maf"
log "		p-value threshold: $p_val"
log "		Genes: $genes"
log "		SNPs: $snps"
log "		Log: $log_file"
log "	Opened $infile"

# Gene selection
if [ "$genes" == "all" ] && [ ! -f "$genes" ]; then
	log "	All genes selected"
	awk -v gidx="$gene_name_idx" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list" \
		|| log "main.sh Error: unable to create gene list" && exit 1
else
	log "	Genes selected from $genes"
	cat "$genes" > "$gene_list"
fi

# SNP selection
if [ "$snps" == "all" ] && [ ! -f "$snps" ]; then
	log "	All SNPs selected"
else
	log "	$snps SNPs selected"
fi

# SNP files
touch snp_count.txt
touch temp_snp_count.txt

# Check if $snps is a file and filter SNPs accordingly
if [ -f "$snps" ]; then
    log "	Filtering SNPs based on $snps"
    if ! awk -v sidx="$snp_id_idx" 'NR==FNR {snps[$1]; next} $sidx in snps' "$snps" "$infile" | sort | uniq -c > temp_snp_count.txt; then
        log "main.sh Error: unable to create filtered snp count file" && exit 1
    fi
else
    if ! awk -v sidx="$snp_id_idx" '{ print $sidx }' "$infile" | sort | uniq -c > temp_snp_count.txt; then
        log "main.sh Error: unable to create snp count file" && exit 1
    fi
fi

if ! awk '{ print $2, $1 }' temp_snp_count.txt > snp_count.txt; then
    log "main.sh Error: unable to create snp_count.txt file" && exit 1
fi

rm temp_snp_count.txt
log "	Created file snp_count.txt"

mkdir -p "$gene_dir"
log "	Created directory $gene_dir"

while IFS= read -r line; do
	log "Working on gene $line..."
	./transform.sh "$line" \
		"$infile" \
		"$gene_dir" \
		"$snp_id_idx" \
		"$chr_idx" \
		"$allele_one_idx" \
		"$allele_two_idx" \
		"$freq_idx" \
		"$gene_name_idx" \
		"$effect_size_idx" \
		"$se_idx" \
		"$p_value_idx" \
		"$snps" \
		"$chr" \
		"$log_file" \
		"$log_dir"
	log "	.ma for $line transformed"
	log "Starting run.sh for $line..."
	./run.sh "$line" \
		"$bfile" \
		"$chr" \
		"$maf" \
		1 \
		"$p_val" \
		"$log_file" \
		"$log_dir"
	log "run.sh for $line finished"
done < "$gene_list" || log "main.sh Error: unable to read gene list" && exit 1
log "main.sh finished"