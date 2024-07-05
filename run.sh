#!/bin/bash

# $1 gene $2 - bfile, $3 - chr, $4 - maf, $5 - run index, $6 p-value threshold
# outputs .cma file
#

gene_name="$1"
bfile="$2"
chr="$3"
maf="$4"
idx="$5"
p_val="$6"
prev_idx="$((idx - 1))"
next_idx="$((idx + 1))"
outfile=$(printf "%s_%s" "$1" "$idx")
infile=$(printf "%s_input.ma" "$1")

# If idx = 1, it means .ma file is used to fetch the lowest p-value
if [ $idx -eq 1 ]
then
	read_file=$infile
	snp_col=1
	p_col=7
else
	read_file=$(printf "%s_%s.cma" "$1" $prev_idx)
	snp_col=2
	p_col=13
fi


top_snp_file=$(printf "%s_%s.snplist" "$1" "$idx")
touch "$top_snp_file"

awk -v col="$p_col" \
	-v id_col="$snp_col" \
	-v thresh="$p_val" \
	'NR == 1 || ($col < thresh && NR > 1 && ($col < min || min == "")) \
	{ min = $col; id = $id_col } END { if (min != "" && min < thresh) \
	print id }' "$read_file" > "$top_snp_file"

# Checks if top snp file is empty
has_snp=$(wc -l < "$top_snp_file")

if [ "$has_snp" -eq 1 ]
then
	gcta64 --bfile "$bfile" --chr "$chr" --maf "$maf" --cojo-file "$infile" \
		--cojo-cond "$top_snp_file" --out "$outfile"
	./run.sh "$gene_name" \
		"$bfile" \
		"$chr" \
		"$maf" \
		"$next_idx" \
		"$p_val"
else
	echo "Total SNPs for $gene_name = $prev_idx"
fi

echo $read_file
