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


infile="$1" # name of the input file
bfile="$2" # name of bfiles
maf="$3"
p_val="$4"

echo "Opened $infile"
touch "$gene_list"


awk -v gidx="$gene_name_idx" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list"

touch snp_count.txt
touch temp_snp_count.txt
awk -v sidx="$snp_id_idx" '{ print $sidx }' "$infile" | sort | uniq -c > \
       temp_snp_count.txt	

awk '{ print $2, $1 }' temp_snp_count.txt > snp_count.txt # SNP: Count file
rm temp_snp_count.txt

mkdir -p "$gene_dir"
while IFS= read -r line; do
	echo "Working on gene $line"
	./transform.sh "$line" \ # gene name
		"$infile" \ # file to be transformed to .ma
		"$gene_dir" \ # name of a directory where data is saved
		"$snp_id_idx" \ # idx of SNP ID in input file
		"$chr_idx" \ # idx of CHR in input file
		"$allele_one_idx" \ # idx of A1 in input file
		"$allele_two_idx" \ # idx of A2 in input file
		"$freq_idx" \ # idx of frequency in input file
		"$gene_name_idx" \ # idx of gene name in input file
		"$effect_size_idx" \ # idx of effect size in input file
		"$se_idx" \ # idx of standard error in input file
		"$p_value_idx" # idx of p-value in input file
	echo ".ma for $line transformed"
	./run.sh "$line" \ # gene name
		"$bfile" \ # bfile
		"$chr" \ # chromosome number
		"$maf" \ # minor allele frequency
		1 \ # first run index 
		"$p_val" # p-value threshold



done < "$gene_list"
echo "Done"
