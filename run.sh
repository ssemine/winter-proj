#!/bin/bash

# $1 gene $2 - bfile, $3 - chr, $4 - maf, $5 .ma file, $6 .snplist file
# $7 run idx
# outputs .cma file
#


idx=$7
prev_idx="$((idx - 1))"
next_idx="$((idx + 1))"

outfile=$(printf "%s_%s" "$1" "$idx")
if [ $idx -eq 1 ]
then
	read_file=$4
	is_cma=0
else
	read_file=$(printf "%s_%s.cma" "$1" $prev_idx)
	is_cma=1
fi

if [ $is_cma -eq 0]
then
	snp_col=1
	p_col=7
else
	snp_col=2
	p_col=13
fi

top_snp_file=$(printf "%s_%s.snplist" "$1" "$idx")
touch "$top_snp_file"

awk -v col="$p_col" -v id_col="$snp_col" -v thresh="$threshold" \
	'NR == 1 || ($col < thresh && NR > 1 && ($col < min || min == "")) \
	{ min = $col; id = $id_col } END { if (min != "" && min < thresh) \
	print id }' "$read_file" > "$top_snp_file"

has_snp=$(wc -l)

if [ $has_snp -eq 1]
then
	gcta64 --bfile $2 --chr $3 --maf $4 --cojo-file $5 \
		--cojo-cond $top_snp_file --out $outfile
	./run.sh $1 $2 $3 $4 $5 $6 $next_idx
else
	echo "Total SNPs for gene = smth"
fi








echo $read_file
