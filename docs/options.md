# User Manual for SLEMM

## Input
### Phenotype file format
- CSV file with a header line
- Required
- The first column must list individual IDs, and each of the following columns lists phenotypes for a trait.
### Genotype file format
- PLINK bed/bim/fam files
- Required
- Family IDs are not used.
### SNP info file format
- CSV file with a header line
- Required
- The first column must list SNP IDs, each of the following columns (if any) lists user-specified prior weights for SNP effect variance.
- Only SNPs available in both the SNP info file and the genotype file are included in the analysis. 
### Covariate file format
- CSV file with a header line
- Optional
- The first column must list individual IDs.

## Options
### REML or LMM
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--phenotype_file` | FILE | Required | Phenotype file |
| `--trait` | STRING` | Required | Must be a column header of the phenotype file and specify which trait to analyze |
| `--error_weight_name` | STRING | Optional | Must be a column header of the phenotype file and specify which column to use for weighting individual error variances |
| `--bfile` | FILE | Required | PLINK bed/bim/fam filename prefix |
| `--snp_info_file` | FILE | Required | SNP info file |
| `--snp_weight_name` |  STRING | Optional | Must be a column header of the SNP info file and specify which column to use for weighting SNP effect variances |
| `--beta_weight_parameters` | STRING | Optional | Specify two parameters (a and b) in the [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) PDF to use scaled beta_pdf(MAF; a, b) as SNP weights  [default is `1,1`] |
| `--covariate_file` | FILE | Optional | Covariate file |
| `--covariate_names` | STRING | Optional | Comma separated list of covariates to include in the analysis |
| `--output_file` | FILE | Required | Output filename or prefix |
| `--window_size` | INT | Optional | Number of SNPs in a window for `--iter_weighting` [default=20] and `--lmm` [default=1000] |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |
| `--subset_size` | INT | Optional | Number of SNPs in a subset [default=1000] |
### SNP filtering
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--min_maf` | FLOAT | Optional | Filter out SNPs with a minor allele frequency below or equal to the provided threshold [default=0] |
| `--min_hwe_pval` | FLOAT | Optional | Filter out SNPs with Hardy-Weinberg equilibrium exact test p-value below the provided threshold [default=0] |
| `--hwe_midp` | FLAG | Optional | Specify the mid-p adjustment in Hardy-Weinberg equilibrium exact tests |
### Stochastic Lanczos REML
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--rel_tol` | FLOAT | Optional | Relative tolerance for the Lanczos decomposition [default=5e-4] |
| `--num_random_vectors` | INT | Optional | Number of random probing vectors [default=30] |
### Predicting genomic breeding values
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--prediction` | FLAG | Required | To predict genomic breeding values |
| `--bfile` | FILE | Required | PLINK bed/bim/fam filename prefix |
| `--snp_estimate_file` | FILE | Required | SNP effect estimate file (e.g., the reml.snp.csv file produced by `--reml` or `--lmm`) |
| `--output` | FILE | Required | Output file where column 1 is individual ID and column 2 is genomic estimated breeding value  |

- This function is similar to --score in PLINK.



## GREML
* --reml
* 

## Genomic Prediction

## GWAS

## Error weight
The argument ```--error_weight_name``` can be added in --reml/--wls/--lmm to weight individual error variance. To enable this, add one additional column in the phenotype CSV file and specify the header name by --error_weight_name. Values are typically 1/r^2-1, where r^2 is the reliability of deregressed EBVs. 
