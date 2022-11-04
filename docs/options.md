## Input
<mark>Missing values in any input file should be left blank. Do not use -9, NA, or NaN.</mark>
### Phenotype file format
- CSV file with a header line
- Required
- The first column must list individual IDs, and each of the following columns lists phenotypes for a trait.
- If phenotypes are deregressed EBVs, an additional column of $(1-r^2)/r^2$ may be available for each trait to represent individual error variance weights, where $r^2$ is the reliability of individual EBVs for that trait. 
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
SLEMM has the following three major functions:
- REML and prediction of SNP effects (`--reml`)
- Mixed-model associations (`--lmm`)
- Prediction of genomic breeding values (`--pred`)

### Common options of REML and LMM
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--phenotype_file` | FILE | Required | Phenotype file |
| `--trait` | STRING` | Required | Must be a column header of the phenotype file and specify which trait to analyze |
| `--error_weight_name` | STRING | Optional | Must be a column header of the phenotype file and specify which column to use as individual error variance weights |
| `--bfile` | FILE PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_info_file` | FILE | Required | SNP info file |
| `--snp_weight_name` |  STRING | Optional | Must be a column header of the SNP info file and specify which column to use as SNP effect variance weights |
| `--beta_weight_parameters` | STRING | Optional | Specify $\alpha$ and $\beta$ as in $f(x;\alpha,\beta)=\text{const}\cdot x^{(\alpha-1)}(1-x)^{(\beta-1)}$ ([beta distribution](https://en.wikipedia.org/wiki/Beta_distribution)) to use $f(\text{MAF};\alpha,\beta)$ as SNP effect variance weights [default is `1,1`] |
| `--covariate_file` | FILE | Optional | Covariate file |
| `--covariate_names` | STRING | Optional | Comma separated list of covariates to include in the analysis |
| `--output_file` | FILE PREFIX | Required | Output filename prefix |
| `--max_heritability` | FLOAT | Optional | Specify a sufficiently large value $h^2_{max}$ so that the search for heritability falls in (1e-4, $h^2_{max}$) [default=0.7] |
| `--lrt` | FLAG | Optional | Flag to perform likelihood-ratio test |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |
| `--window_size` | INT | Optional | Number of SNPs in a window for `--iter_weighting` [default=20] and `--lmm` [default=1000] |
| `--min_maf` | FLOAT | Optional | Filter out SNPs with an MAF below or equal to the provided threshold [default=0] |
| `--min_hwe_pval` | FLOAT | Optional | Filter out SNPs with Hardy-Weinberg equilibrium exact test p-value below the provided threshold [default=0] |
| `--hwe_midp` | FLAG | Optional | Specify the mid-p adjustment in Hardy-Weinberg equilibrium exact tests |
| `--rel_tol` | FLOAT | Optional | Relative tolerance for the Lanczos decomposition [default=5e-4] |
| `--num_random_vectors` | INT | Optional | Number of random probing vectors [default=30] |
| `--seed` | INT | Optional | Random seed [default=0] |
| `--subset_size` | INT | Optional | Number of SNPs in a subset [default=1000] |

### Options specific to REML
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--reml` | FLAG | Required | To run GREML |
| `--iter_weighting` | FLAG | Optional | Flag to run iterative SNP weighting |

When `--iter_weighting` is set, SLEMM has two more options to quickly identify significant, independent SNPs to be fitted as fixed effects. No SNPs can be fitted as fixed by default, because the default chi-square threshold is set to 1e6.

| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--sig_chisq` | FLOAT | Optional | Specify a chi-square threshold of significance for approximate single-SNP associations [default=1e6] |
| `--indep_r2` | FLOAT | Optional | Specify an r2 threshold for identifying independent SNPs from significant ones [default=0.5] |

### Options specific to LMM
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--lmm` | FLAG | Required | To run step 1 for mixed-model associations |
| `--num_qf_markers` | INT | Optional | Specify the number of SNPs for computing quadratic-form statistics [default=30] |

### Options for prediction of genomic breeding values
- This function is similar to --score of PLINK.

| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--prediction` | FLAG | Required | To predict genomic breeding values |
| `--bfile` | FILE PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_estimate_file` | FILE | Required | SNP effect estimate file (e.g., the reml.snp.csv file produced by `--reml` or `--lmm`) |
| `--output` | FILE | Required | Output file where column 1 is individual ID and column 2 is genomic estimated breeding value  |

## Genomic prediction

## Genome-wide association study

