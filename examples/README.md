# Data
A test dataset is available in the **data** folder. 
- **PLINK files**: Genotypes are simulated using [genosim](https://www.ars.usda.gov/northeast-area/beltsville-md-barc/beltsville-agricultural-research-center/agil/aip/software/genosim/).
- **QTL effects** (true_effects.csv): Effects are simulated by sampling from the standard normal distribution.
- **Phenotype files** (.slemm.csv and .bolt.txt): Phenotypes are simulated by summing QTL effects and normal errors, with a heritability of 0.3.

# SLEMM
> [!NOTE]
> Standalone binary executables are available in the [Releases](https://github.com/jiang18/slemm/releases/latest).  
## SNP info file
```console
# Ensure the "data" folder is in your current working directory.
mkdir slemm
cd slemm
perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < ../data/10k.bim > snp_info.csv
```
- The SNP info file is a CSV file with a header line. The file can contain just a single column of SNP IDs.
- The SNP info file specifies which SNPs are included in the random-effects term (or genomic relationship matrix) of the mixed model.
> [!WARNING]
> - **Avoid** listing all sequence variants in the SNP info file when sequence data is available.    
> - **Instead**, list medium-density (e.g., 50k) chip SNPs or LD-pruned variants.
## REML
```console
slemm --reml --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
```
## Window-based SNP weighting
```console
slemm --reml --iter_weighting --window_size 20 --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k.ww --num_threads 10
```
## Prediction of genomic breeding values
```console
# SNP effect estimate file from REML: 10k.reml.snp.csv
slemm --pred --bfile ../data/10k --snp_estimate 10k.reml.snp.csv --out 10k.gebv.csv

# SNP effect estimate file from window-based SNP weighting: 10k.ww.reml.snp.csv
slemm --pred --bfile ../data/10k --snp_estimate 10k.ww.reml.snp.csv --out 10k.ww.gebv.csv
```
## LMM and GWA
```console
# Step 1
slemm --lmm --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
# Step 2
OMP_NUM_THREADS=1 slemm_gamma --pfile ../data/10k --slemm 10k --out 10k.gamma.txt
OMP_NUM_THREADS=10 slemm_gwa --pfile ../data/10k --slemm 10k --out 10k.chr1.txt --chr 1
```
- `slemm --lmm` fits the null linear mixed model for genome-wide associations.
- `slemm_gamma` and `slemm_gwa` use the `--slemm` option to take the output from `slemm --lmm`. 
- `slemm_gamma` computes GRAMMAR-Gamma association statistics for individual SNPs.
- `slemm_gwa` computes single-SNP association statistics that closely approximate those from EMMAX or GCTA-MLMA. 

> [!NOTE]
> Association tests with `slemm_gwa` should be done for each chromosome separately.
> The example below shows how to run tests for chromosomes 1-20 and combine results:
```bash
export OMP_NUM_THREADS=10
for i in `seq 1 20`; do slemm_gwa --pfile ../data/10k --slemm 10k --out 10k.chr$i.txt --chr $i; done

mv 10k.chr1.txt 10k.chrAll.txt
for i in `seq 2 20`; do tail -n +2 10k.chr$i.txt >> 10k.chrAll.txt; rm 10k.chr$i.txt; done
```

# [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
```console
# Ensure the "data" folder is in your current working directory.
mkdir bolt
cd bolt
~/software/BOLT-LMM_v2.3.4/bolt \
  --lmmInfOnly --bed=../data/10k.bed --bim=../data/10k.bim --fam=../data/10k.fam \
  --phenoFile=../data/10k.bolt.txt --phenoCol=QT \
  --numThreads=10 --Nautosomes=30 --LDscoresUseChip \
  --verboseStats --statsFile=10k.lmm.txt
```

# [PLINK](https://www.cog-genomics.org/plink/1.9/)
```console
# Ensure the "data" folder is in your current working directory.
mkdir plink
cd plink
plink --assoc --bfile ../data/10k --allow-no-sex --chr-set 30 --pheno ../data/10k.bolt.txt --pheno-name QT --out 10k
```
