#!/bin/bash

# transform.sh
# ------------
# Transforms GWAS summary statistics file to .ma file format or GCTA's .cma.cojo file to .ma file format.
# File indices are defined in definitions/file_indices.sh.
# -------------------------------------------------------------------------------------------------------
# Usage: ./transform.sh gene_name infile gene_dir snps chr_num log_file log_dir bfile file_type path_to_defitions [ma_file idx]
# -----------------------------------------------------------------------------------------------------------------------------

# Local variables are set
gene_name="$1"
infile="$2"
gene_dir="$3"
chr_num="$4"
log_file="$5"
log_dir="$6"
bfile="$7"
file_type="$8"
path_to_definitions="$9"

# Variables are sourced 
# NOTE: Since transform.sh is called from gene_dir, the path to definitions is relative to gene_dir
source "$path_to_definitions/constants.sh"
source "$path_to_definitions/file_indices.sh"
source "$path_to_definitions/log_messages.sh"
source "$path_to_definitions/functions.sh"

# Checks if transformed file is .ma or .cma.cojo
if [[ "$file_type" != "$INPUT_IDENTIFIER" && "$file_type" != "$CMA_IDENTIFIER" ]]; then
    log_genes "$ERROR_FILE_TYPE"
    exit 1
fi

# Set variables based on file type
if [ "$file_type" = "$CMA_IDENTIFIER" ]; then
    ma_file="${10}"
    idx="${11}"
    name="$(printf "$MA_FILE_NAME_IDX" "$gene_name" "$idx")"
    name_tmp="$(printf "$MA_FILE_NAME_TMP_IDX" "$gene_name" "$idx")"
else
    snps="${10}"
    name="$(printf "$MA_FILE_NAME" "$gene_name")"
    name_tmp="$(printf "$MA_FILE_NAME_TMP" "$gene_name")"
fi

# Set columns for .ma output file
columns="$MA_FILE_COLUMNS"
sample_size=$(wc -l < "$bfile.fam")

# Checks if output files were successfully created
if ! touch "$gene_dir/$name"; then
    log_genes "$ERROR_TOUCH_FAILED $gene_dir/$name"
    exit 1
fi
if ! touch "$gene_dir/$name_tmp"; then
    log_genes "$ERROR_TOUCH_FAILED $gene_dir/$name_tmp"
    exit 1
fi

log_genes "$LOG_WRITING_DATA $gene_dir/$name"

# Transform input file (infile/cma) into .ma file format
if [ "$file_type" = "$INPUT_IDENTIFIER" ]; then
    # If SNPs file is provided
    log "INPUT"
    if [ -f "$snps" ]; then
        log_genes "$LOG_USING_SNP_FILTER $snps"
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
            -v sample_size="$sample_size" \
            'BEGIN {
                while ((getline < snps) > 0) {
                    snp_list[$1]
                }
            }
            {
                if ($gidx == gene && ($col1 in snp_list) && $cidx == chr_num) {
                    print $col1, $col2, $col3, $col4, $col5, $col6, $col7, sample_size >> (gene_dir "/" name)
                }
            }' "$infile" \
            || { log_genes "$ERROR_AWK_WRITE $gene_dir/$name"; exit 1; }
    else
        # If no SNPs file is provided
        log_genes "$LOG_NO_SNP_FILTER"
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
            -v chr_num="$chr_num" \
            -v sample_size="$sample_size" \
            '{
                if (($gidx == gene) && ($cidx == chr_num)) {
                    print $col1, $col2, $col3, $col4, $col5, $col6, $col7, sample_size >> (gene_dir "/" name)
                }
            }' "$infile" \
            || { log_genes "$ERROR_AWK_WRITE $gene_dir/$name"; exit 1; }
    fi
else
    # If .cma.cojo file is provided
    log "CMA"
    awk -F' ' -v col1="$CMA_SNP_ID_IDX" \
    -v col2="$MA_A1_IDX" \
    -v col3="$MA_A2_IDX" \
    -v col4="$CMA_FREQ_IDX" \
    -v col5="$CMA_EFFECT_SIZE_IDX" \
    -v col6="$CMA_SE_IDX" \
    -v col7="$CMA_P_VALUE_IDX" \
    -v gene="$gene_name" \
    -v gene_dir="$gene_dir" \
    -v name="$name" \
    -v ma_snp_idx="$MA_SNP_ID_IDX" \
    -v sample_size="$sample_size" \
    'FNR == 1 { next }
    FNR==NR {
        ma_col2[FNR] = $col2
        ma_col3[FNR] = $col3
        next
    }
    {
        print $col1, ma_col2[FNR], ma_col3[FNR], $col4, $col5, $col6, $col7, sample_size >> (gene_dir "/" name)
    }' "$ma_file" "$infile" \
    || { log_genes "$ERROR_AWK_WRITE $gene_dir/$name"; exit 1; }
fi

log_genes "$LOG_DATA_WRITTEN $gene_dir/$name"

# Final output file transformation
sort -k "$MA_SNP_ID_IDX" "$gene_dir/$name" -o "$gene_dir/$name" \
    || { log_genes "$ERROR_SORT $gene_dir/$name"; exit 1; }
cat "$gene_dir/$name" > "$gene_dir/$name_tmp"
echo "$columns" > "$gene_dir/$name"
cat "$gene_dir/$name_tmp" >> "$gene_dir/$name"
rm "$gene_dir/$name_tmp"

log_genes "$LOG_COLUMNS_ADDED $gene_dir/$name"
log_genes "$LOG_END_MESSAGE $gene_dir/$name"