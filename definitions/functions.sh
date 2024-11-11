# functions.sh
# ------------
# Stores function declarations
# ----------------------------




# General log function
log() {
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo -e "$file_name:$line_number - $message" >> "$log_dir/$log_file"
}

# Log function with gene
log_genes() {
    local gene="$gene_name"
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo -e "$file_name:$line_number - $gene: $message" >> "$log_dir/$log_file"
}

# Log function adding spaces
log_lines() {
    local num_lines="$1"
    for ((i = 0; i < num_lines; i++)); do
        echo "" >> "$log_dir/$log_file"
    done
}

# Log function writing to the summary file
summary_log() {
	local message="$1"
	echo -e "$message" >> "$log_dir/$summary_file"
}