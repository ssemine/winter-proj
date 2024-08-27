#!/bin/bash

dir_name="$1"
outfile="$2"

result_files=$(find "$dir_name"* -name "results")
awk 'FNR==1 && NR!=1 {next;}{print}' $result_files > "$outfile"