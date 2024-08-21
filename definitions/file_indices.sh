# file_indices.sh
# ---------------
# Stores the column indices for the input files: GWAS summary statistics, .ma files, and .cma files.
# --------------------------------------------------------------------------------------------------






# Input file indices (GWAS summary statistics)
export INPUT_SNP_ID_IDX=1
export INPUT_CHR_IDX=2
export INPUT_POS_IDX=3
export INPUT_ALLELE_ONE_IDX=4
export INPUT_ALLELE_TWO_IDX=5
export INPUT_FREQ_IDX=6
export INPUT_GENE_NAME_IDX=7
export INPUT_EFFECT_SIZE_IDX=12
export INPUT_SE_IDX=13
export INPUT_P_VALUE_IDX=14

# .ma file indices
export MA_SNP_ID_IDX=1
export MA_A1_IDX=2
export MA_A2_IDX=3
export MA_FREQ_IDX=4
export MA_EFFECT_SIZE_IDX=5
export MA_SE_IDX=6
export MA_P_VALUE_IDX=7
export MA_SAMPLE_SIZE_IDX=8

# .cma file indices
export CMA_CHR_IDX=1
export CMA_SNP_ID_IDX=2
export CMA_POS_IDX=3
export CMA_ALLELE_EFFECT_IDX=4
export CMA_FREQ_IDX=5
export CMA_ORIGINAL_EFFECT_SIZE_IDX=6
export CMA_ORIGIAL_SE_IDX=7
export CMA_ORIGINAL_P_VALUE_IDX=8
export CMA_ESTIMATED_EFFECTIVE_SAMPLE_SIZE_IDX=9
export CMA_ALLELE_EFFECT_FREQ_IDX=10
export CMA_EFFECT_SIZE_IDX=11
export CMA_SE_IDX=12
export CMA_P_VALUE_IDX=13