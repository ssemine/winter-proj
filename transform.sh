#!/bin/bash

# This script transforms input file to .ma file. 

# Command line arguments indices
gene_name="$1"
infile="$2"
gene_dir="$3"
snp_id_idx="$4"
chr_idx="$5"
allele_one_idx="$6"
allele_two_idx="$7"
freq_idx="$8"
gene_name_idx="$9"
effect_size_idx="$10"
se_idx="$11"
p_value_idx="$14"

# Variables declarations
name="$gene_name.ma"
name_final="${gene_name}_input.ma"
columns="SNP A1 A2 freq b se p N"

# Creates temporary and final .ma files
touch "$gene_dir/$name"
touch "$gene_dir/$name_final"

# Writes required information from input column to .ma file if gene name matches
awk -v "col1=$snp_id_idx" \ 
	-v "col2=$allele_one_idx" \
	-v "col3=$allele_two_idx" \
	-v "col4=$freq_idx" \
	-v "col5=$effect_size_idx" \
	-v "col6=$se_idx" \
	-v "col7=$p_value_idx" \
	-v "gidx=$gene_name_idx" \ # index of gene name
	-v "gene_name=$gene_name"\ # gene name to be matched with
	'{
		if ($gidx == gene_name) {
			print $col1, $col2, $col3, $col4, $col5, $col6, $col7
		}
	}' "$infile" > "$gene_dir/$name"

# Sorts the temporary file by SNP
sort -k 1 "$gene_dir/$name"

# Adds SNP count to N column, writes data to a new file
awk 'NR==FNR {data[$1] = $0; next} $1 in data {print data[$1], $2}' \
	"$gene_dir/$name" snp_count.txt > "$gene_dir/$name_final"

# Removes temporary .ma file
rm "$gene_dir/$name"

# Adds column headers to .ma file
awk -v "cols=$columns" 'BEGIN{print cols}1' "$gene_dir/$name_final" > temp
mv temp "$gene_dir/$name_final"

