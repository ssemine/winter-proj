#!/bin/bash

source definitions/file_indices.sh

gene_name="$1"
infile="$2"
gene_dir="$3"
snps="$4"
chr_num="$5"
log_file="$6"
log_dir="$7"
bfile="$8"
file_type="$9"


log() {
    local gene="$gene_name"
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo "$file_name:$line_number - $gene: $message" >> "$log_dir/$log_file"
}
log_lines() {
    local num_lines="$1"
    for ((i = 0; i < num_lines; i++)); do
        echo "" >> "$log_dir/$log_file"
    done
}
if [[ "$file_type" != "input" && "$file_type" != "cma" ]]; then
    log "Error: file_type must be either 'input' or 'cma'"
    exit 1
fi

if [[ "$file_type" = "cma" ]]; then
    ma_file="${10}"
    run_idx="${11}"
    name=$(printf "%s_%s.ma" "$gene_name" "$idx")
    name_final=$(printf "%s_%s_input.ma" "$gene_name" "$idx")
else
    name=$(printf "%s.ma" "$gene_name")
    name_final=$(printf "%s_input.ma" "$gene_name")
fi

columns="SNP A1 A2 freq b se p N"
sample_size=$(wc -l < "$bfile.fam")

# Creates a two files to put data into
if ! touch "$gene_dir/$name"; then
    log "Error: touch could not create $gene_dir/$name"
    exit 1
fi

if ! touch "$gene_dir/$name_final"; then
    log "Error: touch could not create $gene_dir/$name_final"
    exit 1
fi

# Writes data to the temporary file
log "Writing data to $gene_dir/$name"
if [[ "$file_type" = "input"]]; then
    if [ -f "$snps" ]; then
        log "Using $snps to filter SNPs"
        awk -F' ' -v col1="$INPUT_SNP_ID_IDX" \
            -v col2="$INPUT_ALLELE_ONE_IDX" \
            -v col3="$INPUT_ALLELE_TWO_IDX" \
            -v col4="$INPUT_FREQ_IDX" \
            -v col5="$INPUT_EFFECT_SIZE_IDX" \
            -v col6="$INPUT_SE_IDX" \
            -v col7="$INPUT_P_VALUE_IDX" \
            -v gidx="$INPUT_GENE_NAME_IDX" \
            -v cidx="$INPUT_CHR_IDX" \
            -v gene="$gene_name" \
            -v gene_dir="$gene_dir" \
            -v name="$name" \
            -v snps="$snps" \
            -v chr_num="$chr_num" \
            'BEGIN {
                while ((getline < snps) > 0) {
                    snp_list[$1]
                }
            }
            {
                if ($gidx == gene && ($col1 in snp_list) && $cidx == chr_num) {
                    print $col1, $col2, $col3, $col4, $col5, $col6, $col7 >> (gene_dir "/" name)
                }
            }' "$infile" \
            || { log "Error: awk could not write data to $gene_dir/$name"; exit 1; }
    else
        log "No SNP filter"
        awk -F' ' -v col1="$INPUT_SNP_ID_IDX" \
            -v col2="$INPUT_ALLELE_ONE_IDX" \
            -v col3="$INPUT_ALLELE_TWO_IDX" \
            -v col4="$INPUT_FREQ_IDX" \
            -v col5="$INPUT_EFFECT_SIZE_IDX" \
            -v col6="$INPUT_SE_IDX" \
            -v col7="$INPUT_P_VALUE_IDX" \
            -v gidx="$INPUT_GENE_NAME_IDX" \
            -v cidx="$INPUT_CHR_IDX" \
            -v gene="$gene_name" \
            -v gene_dir="$gene_dir" \
            -v name="$name" \
            -v name_chr="${gene_name}_chr.txt" \
            '{
                if ($gidx == gene) {
                    print $col1, $col2, $col3, $col4, $col5, $col6, $col7 >> (gene_dir "/" name)
                    print gene_name, $cidx >> (gene_dir "/" name_chr)
                }
            }' "$infile" \
            || { log "Error: awk could not write data to $gene_dir/$name"; exit 1; }
    fi
else # need to add that it reads allele_one_idx and two from .ma file but the rest from cma
    awk -F' ' -v col1="$CMA_SNP_ID_IDX" \
        -v col2="$MA_A1_IDX" \
        -v col3="$MA_A2_IDX" \
        -v col4="$CMA_FREQ_IDX" \
        -v col5="$CMA_EFFECT_SIZE_IDX" \
        -v col6="$CMA_SE_IDX" \
        -v col7="$CMA_P_VALUE_IDX" \
        -v col8="$sample_size" \
        -v gene="$gene_name" \
        -v gene_dir="$gene_dir" \
        -v name="$name" \
        'FNR==NR {
            ma_col2[FNR] = $col2
            ma_col3[FNR] = $col3
            next
        }
        {
            if ($gidx == gene) {
                print $col1, ma_col2[FNR], ma_col3[FNR], $col4, $col5, $col6, $col7, $col8 >> (gene_dir "/" name)
            }
        }' "$ma_file" "$infile" \
        || { log "Error: awk could not write data to $gene_dir/$name"; exit 1; }
fi

log "Data written to $gene_dir/$name"

# Sorts the temporary file by SNP
sort -k 1 "$gene_dir/$name" -o "$gene_dir/$name" \
    || { log "Error: sort could not sort $gene_dir/$name"; exit 1; }


# Adds the sample_size number to each row, writes data to a new file
awk -v sample_size="$sample_size" '{print $0, sample_size}' "$gene_dir/$name" > "$gene_dir/$name_final" \
    || { log "Error: awk could not write data to $gene_dir/$name_final"; exit 1; }
log "Sample size $sample_size added to $gene_dir/$name_final"


# Removes temporary .ma file
rm "$gene_dir/$name"

# Adds column headers to .ma file
awk -v "cols=$columns" 'BEGIN{print cols}1' "$gene_dir/$name_final" > temp \
	|| { log "transform.sh Error: awk could not add columns to $gene_dir/$name_final"; exit 1; }
mv temp "$gene_dir/$name_final"
log "Columns added to $gene_dir/$name_final"
mv $gene_dir/$name_final $gene_dir/$name