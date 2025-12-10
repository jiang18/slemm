## Table of contents
- [Input](#input)
- [Options](#options)
    - [Common options of `--reml` and `--lmm`](#common-options-of---reml-and---lmm)
    - [Options specific to `--reml`](#options-specific-to---reml)
    - [Options specific to `--lmm`](#options-specific-to---lmm)
    - [Options for prediction of genomic breeding values](#options-for-prediction-of-genomic-breeding-values)
    - [Options of `slemm_gamma` and `slemm_gwa`](#options-of-slemm_gamma-and-slemm_gwa)
- [Setting options for genomic predictions](#setting-options-for-genomic-predictions)
- [Setting options for GWAS](#setting-options-for-gwas)

## Input
> [!NOTE]
> Missing values in any input CSV file need to be left empty. Do not use space, -9, NA, or NaN.

### Phenotype file format
- CSV file with a header line
- Required
- The first column must be individual IDs, and each subsequent column lists phenotypic values for a trait.
- If phenotypes are deregressed EBVs, an additional column of $(1-r^2)/r^2$ may be present for each trait to represent individual error variance weights, where $r^2$ is the reliability of individual deregressed EBVs for that trait. Refer to `--error_weight_name`.
### Genotype file format
- PLINK bed/bim/fam files
- Required
- Missing genotypes need to be fully filled. Refer to [filling_bed.R](../scripts/filling_bed.R) for a quick how-to.  
- Family IDs are ignored.
### SNP info file format
- CSV file with a header line
- Required
- The first column must list SNP IDs, and each subsequent column (if any) lists user-specified prior weights for SNP effect variance.
- The SNP info file specifies SNPs modeled in the random-effects term (or genomic relationship matrix) of the mixed model.
- Only SNPs present in both the SNP info file and the genotype files are included in the analysis.
### Covariate file format
- CSV file with a header line
- Optional
- The first column must list individual IDs.
- Intercept handling:
   - If `--covariate_names` is not specified, SLEMM automatically includes an intercept term.
   - If `--covariate_names` is specified, SLEMM does not automatically add an intercept. In this case, users must include a column of 1's in the covariate file and specify it in `--covariate_names` if an intercept is desired.

## Options
`slemm` has the following three major functions:
| Function | Flag |
|----------|----------|
| REML and prediction of SNP effects | `--reml` |
| mixed-model GWAS step 1 | `--lmm` |
| prediction of genomic breeding values | `--pred` |

> [!NOTE]
> Option names may be abbreviated if the abbreviation is unique or is an exact match for some defined option; e.g., `--phenotype` works the same as `--phenotype_file`.

### Common options of `--reml` and `--lmm`
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--phenotype_file` | FILE | Required | Phenotype file |
| `--trait` | STRING | Required | Must be a column header in the phenotype file and specify which trait to analyze |
| `--error_weight_name` | STRING | Optional | Must be a column header in the phenotype file and specify which column contains individual error variance weights |
| `--bfile` | FILE PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_info_file` | FILE | Required | SNP info file listing SNPs modeled in the random-effects term of the mixed model |
| `--snp_weight_name` |  STRING | Optional | Must be a column header in the SNP info file and specify which column contains the weights for individual SNPs' contributions to heritability |
| `--beta_weight_parameters` | STRING | Optional | Specify $\alpha$ and $\beta$ as in $f(x;\alpha,\beta)=\text{const}\cdot x^{(\alpha-1)}(1-x)^{(\beta-1)}$ ([beta distribution](https://en.wikipedia.org/wiki/Beta_distribution)) to use $f(\text{MAF};\alpha,\beta)$ as the weights for individual SNPs' contributions to heritability [default is `1,1`] |
| `--covariate_file` | FILE | Optional | Covariate file |
| `--covariate_names` | STRING | Optional | Comma-separated list of covariates to include in the analysis |
| `--output_file` | FILE PREFIX | Required | Output filename prefix |
| `--max_heritability` | FLOAT | Optional | Specify a sufficiently large value $h^2_{max}$ so that the search for heritability falls in (1e-4, $h^2_{max}$) [default=0.7] |
| `--num_threads` | INT | Optional | Number of computational threads to use [default=1] |
| `--window_size` | INT | Optional | Number of SNPs in a window for `--iter_weighting` [default=20] and `--lmm` [default=1000] |
| `--min_maf` | FLOAT | Optional | Filter out SNPs with an MAF less than or equal to the provided threshold [default=0] |
| `--min_hwe_pval` | FLOAT | Optional | Filter out SNPs with Hardy-Weinberg equilibrium exact test p-value below the provided threshold [default=0] |
| `--hwe_midp` | FLAG | Optional | Specify the mid-p adjustment in Hardy-Weinberg equilibrium exact tests |
| `--rel_tol` | FLOAT | Optional | Relative tolerance for the Lanczos decomposition [default=5e-4] |
| `--num_random_vectors` | INT | Optional | Number of random probing vectors [default=30] |
| `--seed` | INT | Optional | Random seed [default=0] |
| `--subset_size` | INT | Optional | SNP subset size for subset-by-subset computations [default=1000] |

> [!NOTE]
> SLEMM weights SNPs in terms of their individual contributions to heritability rather than their squared allele substitution effects.

> [!NOTE]
> `--beta_weight_parameters` will overwrite the SNP weights of `--snp_weight_name` if both options are set.

### Options specific to `--reml`
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--reml` | FLAG | Required | To run GREML |
| `--lrt` | FLAG | Optional | Flag to perform a likelihood-ratio test and to produce standard errors for fixed effects |
| `--iter_weighting` | FLAG | Optional | Flag to run iterative SNP weighting |

When `--iter_weighting` is set, SLEMM has two more options to quickly identify independent significant SNPs to be fitted as fixed effects.
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--sig_chisq` | FLOAT | Optional | Specify a chi-square threshold of significance for approximate single-SNP associations [default=1e6] |
| `--indep_r2` | FLOAT | Optional | Specify an r2 threshold for identifying independent significant SNPs like `plink --clump-r2` [default=0.5] |
> [!NOTE]
> No SNPs can be fitted as fixed by default, because any SNP's chi-square test statistic should in practice be below the default threshold (1e6).

### Options specific to `--lmm`
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--lmm` | FLAG | Required | To run mixed-model GWAS step 1 |
| `--num_qf_markers` | INT | Optional | Specify the number of SNPs for computing quadratic-form statistics [default=30] |

### Options for prediction of genomic breeding values
This function is similar to `--score` of PLINK.
| Option | Argument | Type | Description |
|-------|-------|-------|--------------|
| `--prediction` | FLAG | Required | To predict genomic breeding values |
| `--bfile` | FILE PREFIX | Required | PLINK bed/bim/fam filename prefix |
| `--snp_estimate_file` | FILE | Required | SNP effect estimate file (i.e., the **.reml.snp.csv** file produced by `--reml` or `--lmm`) |
| `--output` | FILE | Required | Output file where column 1 is individual ID and column 2 is genomic estimated breeding value (GEBV)  |

> [!NOTE]
> The estimated SNP effects correspond to A2 alleles in **\*.reml.snp.csv**. GEBV is equal to the sum of SNP effect times A2 allele count across all SNPs. The `slemm --pred` routine can recognize if A1 is swapped with A2 for any SNPs in the prediction population's bim file relative to the training's.

### Options of `slemm_gamma` and `slemm_gwa`
`slemm_gamma` and `slemm_gwa` perform GWAS step 2 after `slemm --lmm` completes step 1. Check options:
```
slemm_gamma --help 
slemm_gwa --help
```
- `slemm_gamma` and `slemm_gwa` use the `--slemm` option to take the output from `slemm --lmm`. 
- `slemm_gamma` computes GRAMMAR-Gamma association statistics for individual SNPs.
- `slemm_gwa` computes single-SNP association statistics that closely approximate those from EMMAX or GCTA-MLMA.
> [!NOTE]
> See [installation guide](./install.md#use-of-slemm_gammapy-and-slemm_gwapy) for using Python scripts instead of standalone executables.

> [!WARNING]
> The genotype files (either bed/bim/fam or pgen/pvar/psam) taken by `slemm_gwa`'s `--pfile` option must contain all the null-model SNPs on the specified chromosome;
> otherwise, an error will occur. Null-model SNPs are shown in `slemm --lmm`'s `.reml.snp.csv` output file.

## Setting options for genomic predictions
- `--reml` is required.
- `--num_threads` should be proporly set on a multi-core computer to speed up computations.
- `--max_heritability` should be set **slightly larger** than the actual (pseudo-)heritability; e.g., `--max_herit 0.4` can be used for dairy milk yield traits because we know they are moderately heritable. A large value works but results in unnecessary computation. Its default value (0.7) is **unnecessarily** large for most traits.
- `--iter_weighting` is particularly needed for traits underlain by large-effect QTLs.
- `--window_size` is only useful when `--iter_weighting` is set. Its default value generally works well.

## Setting options for GWAS
- `--lmm` is required.
- `--num_threads` should be proporly set on a multi-core computer to speed up computations.
- `--max_heritability` should be set **slightly larger** than the actual (pseudo-)heritability; e.g., `--max_herit 0.4` can be used for dairy milk yield traits because we know they are moderately heritable. A large value works but results in unnecessary computation. Its default value (0.7) is **unnecessarily** large for most traits.
- `--window_size` affects the computation of mixed-model associations in `slemm_gwa`. The larger `--window_size`, the more accurate approximation of mixed-model associations but the more intensive computations. Its default value generally works well for livestock data.
> [!NOTE]
> The behavior of `--window_size` in `--reml` differs from that in `--lmm`.
