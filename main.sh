#!/bin/bash

# TODO: Add creating a file in main.sh that would be a list GENE: CHR
# Fix issues with format strings
# Testing:
# transform.sh -> works
# main.sh, run.sh -> need to test
# Add parallel functionality


# Arguments:
#
# $1 infile, $2 bfile, $3 maf, $4 p-value thresh

# CONSTANTS
snp_id_idx=1
chr_idx=2
allele_one_idx=4
allele_two_idx=5
freq_idx=6
gene_name_idx=7
p_value_idx=14
gene_list="genelist.txt"
gene_dir="cojo_files"


infile="$1"
bfile="$2"
maf="$3"
p_val="$4"

echo "Opened $infile"
touch "$gene_list"


awk -v gidx="$gene_name_idx" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list"

touch snp_count.txt
touch temp_snp_count.txt
awk -v sidx="$snp_id_idx" '{ print $sidx }' "$infile" | sort | uniq -c > \
       temp_snp_count.txt	

awk '{ print $2, $1 }' temp_snp_count.txt > snp_count.txt
rm temp_snp_count.txt

mkdir -p "$gene_dir"
while IFS= read -r line; do
	echo "Working on gene $line"
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
	./run.sh "$line" \
		"$bfile" \
		"$chr" \
		"$maf" \
		1 \
		"$p_val"
done < "$gene_list"
echo "Done"
