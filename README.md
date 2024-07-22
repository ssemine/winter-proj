### Usage

Run 

`./main.sh --infile input_file --bfile bfile --maf maf --p_val p_val --chr chr_num 
[ --genes gene_set | --gene_dir gene_directory | --log log_file | --snps snp_set]`

### Required files:
gcta64 binary in this directory

Input file

Bfiles (.bed, .bim, .fam)


### Arguments:

  Required:
  
    --infile - input file
    --bfile - bed file
    --maf - minor allele frequency
    --p_val - p-value threshold
    --chr - chromosome number
    
  Optional:
  
    --genes - file with subset of genes to run GCTA-COJO for (Row format)
    --gene_dir - directory name where gene .ma files will be saved
    --log - log file name
    --snps - file with subset of SNPs to run GCTA-COJO for (Row format)

### Files

  **main.sh** - main file which should be executed by the user
  
  **transform.sh** - transforms infile (input file) into a .ma file for each gene
  
  **run.sh** - runs GCTA-COJO iteratively for each SNP


