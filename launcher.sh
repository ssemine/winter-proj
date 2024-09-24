#!/bin/bash

infile="$1"
bfiles_dir="$2"
maf="$3"
p_val="$4"
exclude_qtl_type="$5"

bfiles=$(ls "$bfiles_dir" | cut -d. -f1 | sort | uniq)

bim_chr_idx=1
bim_snp_idx=2

original_path=$(pwd)

for bfile in $bfiles; do
    chr=$(echo "$bfile" | grep -oP '(?<=chrom)\d+')
    new_bfiles_dir="$bfiles_dir"
    if [[ "$chr" -ge 23 && "$chr" -le 29 ]]; then
        echo "Changed from $chr"
        bim_chr=$(echo "$chr - 20" | bc)
        echo "to $bim_chr"
	cd "$bfiles_dir"
	cd ..
	new_bfiles_dir=$(pwd)
	if [[ -d "$new_bfiles_dir/$bfile" ]]; then
		rm -rf "$new_bfiles_dir/$bfile"
		echo "Removed $new_bfiles_dir/$bfile"
	fi
	cd "$original_path"
        mkdir "$new_bfiles_dir/$bfile"
        cp "$bfiles_dir/$bfile"* "$new_bfiles_dir/$bfile/."
        new_bfiles_dir="$new_bfiles_dir/$bfile"
	mv "$new_bfiles_dir/$bfile.bim" "$new_bfiles_dir/$bfile.bim.tmp"
        awk -v chr_idx="$bim_chr_idx" \
            -v snp_idx="$bim_snp_idx" \
            -v chr="$bim_chr" \
            '{
                split($snp_idx, snp, ":");
		$(chr_idx) = chr;
		$(snp_idx) = chr ":" snp[2];
                print $0;
            }' "$new_bfiles_dir/$bfile.bim.tmp" > "$new_bfiles_dir/$bfile.bim"
    fi
    sbatch slurm.sh "$infile" \
        "$new_bfiles_dir/$bfile" \
        "$maf" \
        "$p_val" \
        "$chr" \
        "$exclude_qtl_type" \
    	"$bfile"
done
