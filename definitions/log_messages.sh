# log_messages.sh
# ---------------
# Stores the log and error messages.
# ----------------------------------





# MAIN.SH MESSAGES
# ----------------

# error messages
export ERROR_INVALID_ARGUMENT="Error: invalid argument:"
export ERROR_CHR_NUM_NOT_PROVIDED="Error: chromosome number not provided"
export ERROR_GENE_LIST="Error: unable to create gene list"
export ERROR_FILTERED_SNP_COUNT_FILE="Error: unable to create filtered SNP count file"
export ERROR_SNP_COUNT_FILE="Error: unable to create SNP count file"
export ERROR_READ_GENE_LIST="Error: unable to read gene list"

# log messages
export LOG_WELCOME_MESSAGE="\n---------------------\n"
LOG_WELCOME_MESSAGE+="GCTA-COJO script\n"
LOG_WELCOME_MESSAGE+="\n"
LOG_WELCOME_MESSAGE+="---------------------\n"
export LOG_PARAMETERS="Parameters selected:\n"
LOG_PARAMETERS+="   Input file: %s\n"
LOG_PARAMETERS+="   BED files: %s\n"
LOG_PARAMETERS+="   Minor allele frequency: %s\n"
LOG_PARAMETERS+="   P-value threshold: %s\n"
LOG_PARAMETERS+="   Chromosome number: %s\n"
LOG_PARAMETERS+="   Genes: %s\n"
LOG_PARAMETERS+="   SNPs: %s\n"
LOG_PARAMETERS+="   Log file: %s\n"
export LOG_CREATED_FILES="Created file(s):"
export LOG_CALLING_TRANSFORM="Calling transform.sh for"
export LOG_CALLING_RUN="Calling run.sh for"
export LOG_MA_TRANSFORMED=".ma file transformed for"
export LOG_RUN_FINISHED="run.sh finished for"
export LOG_END_MESSAGE="End of script"





# TRANSFORM.SH MESSAGES
# ---------------------

# error messages
export ERROR_FILE_TYPE="Error: file_type must be either 'input' or 'cma'"
export ERROR_TOUCH_FAILED="Error: touch could not create"
export ERROR_AWK_WRITE="Error: awk could not write data to"
export ERROR_SORT="Error: sort failed for"
export ERROR_AWK_COLS="Error: awk could not add columns to"

# log messages
export LOG_STARTING_TRANSFORM="Starting transform.sh for"
export LOG_WRITING_DATA="Writing data to"
export LOG_USING_SNP_FILTER="Using SNP filter:"
export LOG_NO_SNP_FILTER="No SNP filter"
export LOG_DATA_WRITTEN="Data written to"
export LOG_SAMPLE_SIZE_ADDED="Sample size %s added to %s"
export LOG_SNP_COUNT_ADDED="SNP count added to"
export LOG_COLUMNS_ADDED="Columns added to"
export LOG_FILE_RENAMED="%s renamed to %s"
export LOG_ENDING_TRANSFORM="Ending transform.sh for"





# RUN.SH MESSAGES
# ---------------

# error messages

# log messages
export LOG_STARTING_RUN="Starting run.sh for"
export LOG_ITERATION_NUM="Iteration number:"
export LOG_READING_FROM="Reading from"
export LOG_TOP_SNP="Top SNP for %s: %s"
export LOG_TOTAL_SNPS="Total SNPs: %s: %s"
export LOG_ENDING_RUN="Ending run.sh for"
