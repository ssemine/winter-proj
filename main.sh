#!/bin/bash

# TODO: Add creating a file in main.sh that would be a list GENE: CHR
# Fix issues with format strings
# Testing:
# transform.sh -> works
# main.sh, run.sh -> need to test
# Add parallel functionality

# Add optional arguments in --format


# Arguments:
#
# $1 infile, $2 bfile, $3 maf, $4 p-value thresh

# main.sh --infile infile --bfile bfile --maf maf --p-value p_val --gene-names gene_names 

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
        --genes)
            genes="$2"
            shift 2
            ;;
        *)
            echo "Error: invalid argument: $1"
            exit 1
            ;;
    esac
done

echo "Running main.sh..."
echo "Parameters selected: "
echo "	Input file: $infile"
echo "	BED files: $bfile"
echo "	Minor allele frequency: $maf"
echo "	p-value threshold: $p_val"

echo "Opened $infile"


if [ "$genes" == "all" ]
then
	echo "All genes selected"
	awk -v gidx="$gene_name_idx" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list"
else
	echo "Genes selected from $genes"
	echo "$genes" > "$gene_list"
fi

touch snp_count.txt
touch temp_snp_count.txt

awk -v sidx="$snp_id_idx" '{ print $sidx }' "$infile" | sort | uniq -c > \
       temp_snp_count.txt	
awk '{ print $2, $1 }' temp_snp_count.txt > snp_count.txt
rm temp_snp_count.txt
echo "Created snp_count.txt"

mkdir -p "$gene_dir"
echo "Created $gene_dir"
while IFS= read -r line; do
	echo "Working on gene $line..."
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
		"$chr_idx"
	echo ".ma for $line transformed"
	chr=$(grep "^$line " "$gene_dir/${line}_chr.txt" | awk '{print $2}')
	echo "Starting run.sh for $line..."
	./run.sh "$line" \
		"$bfile" \
		"$chr" \
		"$maf" \
		1 \
		"$p_val"
	echo "run.sh for $line finished"
done < "$gene_list"
echo "Done"
