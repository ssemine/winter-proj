### Usage

Run 

`./main.sh [ required args ] [ optional args ]`

Results will be saved in `results` file, with the following headers:

`"SNP CHR BP START GENE STRAND A1 A2 FREQ B B_C SE SE_C P P_C P_T N TYPE ROUND"`

### Required files:
gcta64 binary in this directory (otherwise specify path in definitions/contants.sh)

Input file (SNPs of interest)

Bfiles (.bed, .bim, .fam)


### Arguments:

  Required:
  
    --infile - input file
    --bfile - bed file
    --chr - chromosome number
    --maf - minor allele frequency
    
  Optional:

    
    --genes - file with subset of genes to run GCTA-COJO for (Row format) (gene_list by default)
    --gene_dir - directory name where gene .ma files will be saved (cojo_files by default)
    --run-dir - directory name where all the files will be saved (run_$date by default)
    --log - log file name ($date.log by default)
    --p_val - p-value threshold (value or "per_chr" or "per_gene")
    --snps - file with subset of SNPs to run GCTA-COJO for (Row format)

### Files

  **main.sh** - main file which should be executed by the user
  
  **transform.sh** - transforms infile (input file) into a .ma file for each gene
  
  **run.sh** - runs GCTA-COJO iteratively for each SNP

  **definitions/constants.sh** - stores constants, such as the default p-value

  **definitions/file_indices.sh** - stores input (GWAS summary statistics), .ma, .cma.cojo column indices

  **definitions/functions.sh** - stores logging functions' declarations

  **definitions/log_messages.sh** - stores log messages

  


