#!/bin/bash

# main.sh --infile infile --bfile bfile --maf maf --p-value p_val --gene gene_names 

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
echo "	Parameters selected: "
echo "		Input file: $infile"
echo "		BED files: $bfile"
echo "		Minor allele frequency: $maf"
echo "		p-value threshold: $p_val"

echo "	Opened $infile"


if [ "$genes" == "all" ]
then
	echo "	All genes selected"
	awk -v gidx="$gene_name_idx" '{ print $gidx }' "$infile" | sort | uniq > "$gene_list" \
		|| echo "main.sh Error: unable to create gene list" && exit 1
else
	echo "	Genes selected from $genes"
	echo "$genes" > "$gene_list"
fi

touch snp_count.txt
touch temp_snp_count.txt

awk -v sidx="$snp_id_idx" '{ print $sidx }' "$infile" | sort | uniq -c > \
       temp_snp_count.txt \
	   || echo "main.sh Error: unable to create snp count file" && exit 1

awk '{ print $2, $1 }' temp_snp_count.txt > snp_count.txt \
	|| echo "main.sh Error: unable to create snp_count.txt file" && exit 1

rm temp_snp_count.txt
echo "	Created file snp_count.txt"

mkdir -p "$gene_dir"
echo "	Created directory $gene_dir"
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
	echo "	.ma for $line transformed"
	chr=$(grep "^$line " "$gene_dir/${line}_chr.txt" | awk '{print $2}') \
		|| echo "main.sh Error: unable to fetch chromosome number" && exit 1
	echo "Starting run.sh for $line..."
	./run.sh "$line" \
		"$bfile" \
		"$chr" \
		"$maf" \
		1 \
		"$p_val"
	echo "run.sh for $line finished"
done < "$gene_list" || echo "main.sh Error: unable to read gene list" && exit 1
echo "main.sh finished"