# User Manual for SLEMM

## Input
### Phenotype file format
- CSV file with a header line
- The phenotype file is required.
- The first column must list individual IDs, and each of the following columns lists phenotypes for a trait.
### Genotype file format
- PLINK bed/bim/fam files
- The genotype file is required.
- Family IDs are not used.
### SNP info file format
- CSV file with a header line
- The SNP info file is required.
- The first column must list SNP IDs.
- Only SNPs available in both the SNP info file and the genotype file are included in analysis. 
### Covariate file format
- CSV file with a header line
- The covariate file is optional.
- The first column must list individual IDs.

## GREML
* --reml
* 

## Genomic Prediction

## GWAS

## Error weight
The argument ```--error_weight_name``` can be added in --reml/--wls/--lmm to weight individual error variance. To enable this, add one additional column in the phenotype CSV file and specify the header name by --error_weight_name. Values are typically 1/r^2-1, where r^2 is the reliability of deregressed EBVs. 
