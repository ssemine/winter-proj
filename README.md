### Usage

Run 

`./main.sh --infile input_file --bfile bfile --chr chr_num --maf maf  
[ --genes gene_set | --gene_dir gene_directory | --log log_file | --p_val p_val | --snps snp_set]`

### Required files:
gcta64 binary in this directory

Input file

Bfiles (.bed, .bim, .fam)


### Arguments:

  Required:
  
    --infile - input file
    --bfile - bed file
    --chr - chromosome number
    --maf - minor allele frequency
    
  Optional:

    
    --genes - file with subset of genes to run GCTA-COJO for (Row format)
    --gene_dir - directory name where gene .ma files will be saved
    --log - log file name
    --p_val - p-value threshold
    --snps - file with subset of SNPs to run GCTA-COJO for (Row format)

### Files

  **main.sh** - main file which should be executed by the user
  
  **transform.sh** - transforms infile (input file) into a .ma file for each gene
  
  **run.sh** - runs GCTA-COJO iteratively for each SNP

  **definitions/constants.sh** - stores constants, such as the default p-value

  **definitions/file_indices.sh** - stores input (GWAS summary statistics), .ma, .cma.cojo column indices

  **definitions/functions.sh** - stores logging functions' declarations

  **definitions/log_messages.sh** - stores log messages


