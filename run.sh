#!/bin/bash

# run.sh
# ------
# Fetches lowest p-value SNP from .ma file and runs GCTA conditional analysis. 
# Calls transform.sh to transform the GCTA's ouput .cma.cojo to .ma file, passes it to run.sh until all SNPs are with p-value < threshold are fetched.
# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Usage: ./run.sh gene_name bfile chr maf idx p_val log_file log_dir gene_dir snp_dir summary_file
# ------------------------------------------------------------------------------------------------


gene_name="$1"
bfile="$2"
chr="$3"
maf="$4"
idx="$5"
p_val="$6"
log_file="$7"
log_dir="$8"
gene_dir="$9"
snp_dir="${10}"
summary_file="${11}"
path_to_definitions="${12}"
results_file="${13}"


source "$path_to_definitions/constants.sh"
source "$path_to_definitions/file_indices.sh"
source "$path_to_definitions/log_messages.sh"

prev_idx="$((idx - 1))"
next_idx="$((idx + 1))"
outfile="$(printf "$GCTA_OUTFILE_NAME" "$gene_name" "$idx")"
infile="$(printf "$MA_FILE_NAME" "$gene_name")"
ma_file_reference="$(printf "$MA_FILE_NAME_REFERENCE" "$gene_name")"
ma_file_reference_tmp="$(printf "$MA_FILE_NAME_REFERENCE_TMP" "$gene_name")"

source "$path_to_definitions/functions.sh"

log "$LOG_STARTING_RUN $gene_name"
log "$LOG_ITERATION_NUM $idx"

if [ $idx -eq 1 ]; then
	read_file="$infile"
else
	read_file="$(printf "$MA_FILE_NAME_IDX" "$gene_name" $prev_idx)"
fi
log "$LOG_READING_FROM $gene_dir/$read_file"

top_snp_file="$(printf "$TOP_SNP_FILE" "$snp_dir" "$gene_name" "$idx")"
touch "$top_snp_file"

awk -v col="$MA_P_VALUE_IDX" \
    -v id_col="$MA_SNP_ID_IDX" \
    -v thresh="$p_val" \
    -v snp_file="$top_snp_file" \
    -v summary_file="$log_dir/$summary_file" \
    'NR > 1 && $col < thresh { 
        if (min == "" || $col < min) { 
            min = $col; 
            id = $id_col 
        } 
    } 
    END { 
        if (min != "" && min < thresh) 
            print id > snp_file
            print id, min >> summary_file
    }' "$gene_dir/$read_file" \
    || { log "$ERROR_AWK_WRITE $top_snp_file"; exit 1; }

has_snp="$(wc -l < "$top_snp_file")"
touch "$MA_TOP_SNP_FILE"
if [[ "$has_snp" =~ ^-?[0-9]+$ ]] && [ "$has_snp" -eq 1 ]; then
    top_snp=$(cat $top_snp_file)
	log "$(printf "$LOG_TOP_SNP" "$gene_name" "$top_snp")"
	"$PATH_TO_GCTA" --bfile "$bfile" \
        --chr "$chr" \
        --maf "$maf" \
        --cojo-file "$gene_dir/$read_file" \
		--cojo-cond "$top_snp_file" \
        --out "$gene_dir/$outfile" \
        || { log "$ERROR_GCTA_FAILED $gene_name"; exit 1; }
    cma_file="$(printf "$TRANSFORM_CMA_FILE_NAME" "$gene_dir" "$outfile")"
    awk -v snp="$MA_SNP_ID_IDX" -v top_snp="$top_snp" '$snp == top_snp' "$gene_dir/$read_file" > "$MA_TOP_SNP_FILE"
    if [[ "$idx" =~ ^-?[0-9]+$ ]] && [ "$idx" -ge 2 ]; then
        touch "$CMA_TOP_SNP_FILE"
        prev_cma_file="$(printf "$TRANSFORM_CMA_FILE_NAME" "$gene_dir" "$(printf "$GCTA_OUTFILE_NAME" "$gene_name" $prev_idx)")"
        awk -v snp="$CMA_SNP_ID_IDX" -v top_snp="$top_snp" '$snp == top_snp' "$prev_cma_file" > "$CMA_TOP_SNP_FILE"
        awk -v snp="$MA_SNP_ID_IDX" \
            -v gene_name="$gene_name" \
            -v allele_one="$MA_A1_IDX" \
            -v allele_two="$MA_A2_IDX" \
            -v freq="$MA_FREQ_IDX" \
            -v effect_size="$MA_EFFECT_SIZE_IDX" \
            -v se="$MA_SE_IDX" \
            -v p_val="$MA_P_VALUE_IDX" \
            -v cma_effect_size_idx="$CMA_EFFECT_SIZE_IDX" \
            -v cma_se_idx="$CMA_SE_IDX" \
            -v cma_p_val_idx="$CMA_P_VALUE_IDX" \
            -v sample_size="$MA_SAMPLE_SIZE_IDX" \
            -v thresh="$p_val" \
            'FNR==NR {
                cma_effect_size=$cma_effect_size_idx
                cma_se=$cma_se_idx
                cma_p_val=$cma_p_val_idx
                next
            }
            { 
                print $snp, gene_name, $allele_one, $allele_two, $freq, $effect_size, $se, $p_val, cma_effect_size, cma_se, cma_p_val, $sample_size, thresh
            }' "$CMA_TOP_SNP_FILE" "$MA_TOP_SNP_FILE" >> "$results_file"
        rm "$CMA_TOP_SNP_FILE"
    fi

    head -n 1 "$cma_file" > "$cma_file.$HEADER_EXTENTION"
    tail -n +2 "$cma_file" | sort -k "$CMA_SNP_ID_IDX" > "$cma_file.$SORTED_EXTENTION"
    cat "$cma_file.$HEADER_EXTENTION" "$cma_file.$SORTED_EXTENTION" > "$cma_file"
    rm "$cma_file.$HEADER_EXTENTION" "$cma_file.$SORTED_EXTENTION"  

    awk -v snp_col="$CMA_SNP_ID_IDX" \
        'NR > 1 {
            print $snp_col
        }' "$cma_file" > "$gene_dir/$SNP_LIST"

    awk -v snp_col="$MA_SNP_ID_IDX" \
        'NR==FNR { 
            snps[$snp_col];
            next 
        } 
        FNR == 1 { 
            print; 
            next
        } 
        $snp_col in snps' "$gene_dir/$SNP_LIST" \
        "$gene_dir/$ma_file_reference" > "$gene_dir/$ma_file_reference.$TMP_EXTENTION"

    cat "$gene_dir/$ma_file_reference.$TMP_EXTENTION" > "$gene_dir/$ma_file_reference"
    rm "$gene_dir/$SNP_LIST" "$gene_dir/$ma_file_reference.$TMP_EXTENTION"

    "$PATH_TO_TRANSFORM_SH" "$gene_name" \
        "$cma_file" \
        "$gene_dir" \
        "$chr" \
        "$log_file" \
        "$log_dir" \
        "$bfile" \
        "$CMA_IDENTIFIER" \
        "$PATH_TO_DEFINITIONS" \
        "$gene_dir/$ma_file_reference" \
        "$idx" \
        || { log "$ERROR_TRANSFORM $gene_name"; exit 1; }
	"$PATH_TO_RUN_SH" "$gene_name" \
		"$bfile" \
		"$chr" \
		"$maf" \
		"$next_idx" \
		"$p_val" \
		"$log_file" \
		"$log_dir" \
		"$gene_dir" \
        "$snp_dir" \
        "$summary_file" \
        "$PATH_TO_DEFINITIONS" \
        "$results_file" \
        || { log "$ERROR_RUN_FAILED $gene_name"; exit 1; }
else
    if [[ "$prev_idx" =~ ^-?[0-9]+$ ]] && [ "$prev_idx" -eq 1 ]; then
        cat "$MA_TOP_SNP_FILE"
        awk -v snp="$MA_SNP_ID_IDX" \
            -v gene_name="$gene_name" \
            -v allele_one="$MA_A1_IDX" \
            -v allele_two="$MA_A2_IDX" \
            -v freq="$MA_FREQ_IDX" \
            -v effect_size="$MA_EFFECT_SIZE_IDX" \
            -v se="$MA_SE_IDX" \
            -v p_val="$MA_P_VALUE_IDX" \
            -v sample_size="$MA_SAMPLE_SIZE_IDX" \
            -v thresh="$p_val" \
            '{ 
                print $snp, gene_name, $allele_one, $allele_two, $freq, $effect_size, $se, $p_val, $effect_size, $se, $p_val, $sample_size, thresh
            }' "$MA_TOP_SNP_FILE" >> "$results_file"
    fi
	log "$(printf "$LOG_TOTAL_SNPS" "$gene_name" "$prev_idx")"
    summary_log "$(printf "$LOG_TOTAL_SNPS" "$gene_name" "$prev_idx")"
fi