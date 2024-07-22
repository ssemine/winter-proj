### Usage

Run 

`./main.sh --infile input_file --bfile bfile --maf maf --p_val p_val --chr chr_num 
[ --genes gene_set | --gene_dir gene_directory | --log log_file | --snps snp_set]`

### Arguments:

  Required:
  
    --infile - input file
    --bfile - bed file
    --maf - minor allele frequency
    --p_val - p-value threshold
    --chr - chromosome number
    
  Optional:
  
    --genes - file with subset of genes to run GCTA-COJO for. (Row format)
    --gene_dir - directory name where gene .ma files will be saved
    --log - log file name
    --snps - file with subset of SNPs to run GCTA-COJO for. (Row format)
    
    
