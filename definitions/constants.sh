# constants.sh
# ------------
# Stores the constants.
# ---------------------





export P_VALUE_THRESHOLD="1e-10"
export P_VALUE_THRESHOLD_PER_GENE="per-gene"
export P_VALUE_THRESHOLD_PER_CHR="per-chr"
export GENE_LIST="gene_list"
export GENE_DIR="cojo_files"
export RUN_DIR="run_$(date '+%Y-%m-%d %H:%M:%S')"
export SNP_DIR="snp_files"
export GENES_ALL="all"
export SNPS_ALL="all"
export LOG_DIR="logs"
export LOG_FILE="$(date '+%Y-%m-%d %H:%M:%S').log"
export SNP_CSV_HEADER="SNP,Count"
export SNP_COUNT_FILE="snp_count.csv"
export SUMMARY_FILE="summary.log"
export MA_FILE_COLUMNS="SNP A1 A2 freq b se p N"
export MA_FILE_NAME="%s.ma"
export MA_FILE_NAME_FINAL="%s_input.ma"
export MA_FILE_NAME_IDX="%s_%s.ma"
export MA_FILE_NAME_FINAL_IDX="%s_%s_input.ma"
export MA_FILE_NAME_REFERENCE="%s_reference.ma"
export MA_FILE_NAME_REFERENCE_TMP="%s_reference.ma.tmp"
export CMA_IDENTIFIER="cma"
export INPUT_IDENTIFIER="input"
export TOP_SNP_FILE="%s/%s_%s.snplist"
export TRANSFORM_CMA_FILE_NAME="%s/%s.cma.cojo"
export GCTA_OUTFILE_NAME="%s_%s"