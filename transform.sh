#!/bin/bash


# Constants
snp_id_idx=1
chr_idx=2
allele_one_idx=4
allele_two_idx=5
freq_idx=6
gene_name_idx=7
effect_size_idx=12
se_idx=13
p_value_idx=14

n_args=3
if [ "$#" -lt "$n_args" ]
then
	echo "Gene/Infile is not supplied"
	exit 1
fi

gene_name="$1"
infile="$2"
gene_dir="$3"
name="$gene_name.ma" # name per gene?
name_final="${gene_name}_input.ma"
columns="SNP A1 A2 freq b se p N"

touch "$gene_dir/$name"
touch "$gene_dir/$name_final"
# echo "$columns" > "$name" # sets headers for file 
awk -v "col1=$snp_id_idx" -v "col2=$allele_one_idx" -v "col3=$allele_two_idx" \
	-v "col4=$freq_idx" -v "col5=$effect_size_idx" -v "col6=$se_idx" \
	-v "col7=$p_value_idx" -v "gidx=$gene_name_idx" -v "gene_name=$gene_name"\
	'{
		if ($gidx == gene_name) {
			print $col1, $col2, $col3, $col4, $col5, $col6, $col7
		}
	}' "$infile" > "$gene_dir/$name"

sort -k 1 "$gene_dir/$name" # Sorted by SNP

awk 'NR==FNR {data[$1] = $0; next} $1 in data {print data[$1], $2}' \
	"$gene_dir/$name" snp_count.txt > "$gene_dir/$name_final"

rm "$gene_dir/$name"

awk -v "cols=$columns" 'BEGIN{print cols}1' "$gene_dir/$name_final" >tmp
mv tmp "$gene_dir/$name_final"

