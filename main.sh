#!/bin/bash



# We need to transform the input file to a .ma file that is used with --cojo-file. <- separate script
# Format of .ma fike:
# SNP A1 A2 freq b se p N
#
#
# we also need to create a temporary file for storing all the gene names (and maybe respective chromosome)
# then with a list of SNPs -> --cojo-cond cond.snplist
# gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-cond cond.snplist --out test_chr1

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
	./transform.sh "$line" "$infile" "$gene_dir"
done < "$gene_list"	


# Now, I shoud extract an SNP with lowest p-value from resulted .ma file. Then, I run GCTA-COJO that returns .cma 
#
# I select a p-value threshold and until there are no SNPs that are below that threshold, I keep doing GCTA-COJO
# for each. one SNP in .snplist file. Would probably need a while loop to iterate over the .cma file, find the lowest
# and perform GCTA-COJO again and again. Maybe think of  a different way to iterate over all genes and all SNPs for 
# each gene. Also, choose tools that scale well, so there is no computational limit to the script. 
#
# 01/07
