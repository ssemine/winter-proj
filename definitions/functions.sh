log() {
    local message="$1"
    local line_number="${BASH_LINENO[0]}"
    local file_name="${BASH_SOURCE[1]}"
    echo "$file_name:$line_number - $message" >> "$log_dir/$log_file"
}
log_genes() {
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
summary_log() {
	local message="$1"
	echo "$message" >> "$log_dir/$summary_file"
}
