# Data
A test dataset is available in the **data** folder. 
- **PLINK files**: Genotypes are simulated using [genosim](https://aipl.arsusda.gov/software/genosim/).
- **QTL effects** (true_effects.csv): Effects are simulated by sampling from the standard normal distribution.
- **Phenotype files** (.slemm.csv and .bolt.csv): Phenotypes are simulated by summing QTL effects and normal errors, with a heritability of 0.3.

# SLEMM
## SNP info file
```console
# Ensure the "data" folder is in your current working directory.
mkdir slemm
cd slemm
perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < ../data/10k.bim > snp_info.csv
```
- The SNP info file is a CSV file with a header line and specifies what SNPs to be included in SLEMM.
## REML
```console
slemm --reml --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
```
## Prediction of genomic breeding values
```console
slemm --pred --bfile ../data/10k --snp_estimate 10k.reml.snp.csv --out 10k.gebv.csv
```
## LMM and GWA
```console
slemm --lmm --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
OMP_NUM_THREADS=1 slemm_gamma.py --pfile ../data/10k --slemm 10k --out 10k.gamma.txt
OMP_NUM_THREADS=10 slemm_gwa.py --pfile ../data/10k --slemm 10k --out 10k.chr1.txt --chr 1
```
- `slemm --lmm` fits linear mixed model for genome-wide associations.
- `slemm_gamma.py` and `slemm_gwa.py` use the option `--slemm` to take the output of `slemm --lmm`. 
- `slemm_gamma.py` computes GRAMMAR-Gamma association statistics for each SNP.
- `slemm_gwa.py` computes single-SNP association statistics that closely approximate those obtained from EMMAX or GCTA-MLMA. 
> [!NOTE]
> Association tests with `slemm_gwa` should be done for each chromosome separately. 
```bash
export OMP_NUM_THREADS=10
for i in `seq 1 30`; do slemm_gwa.py --pfile ../data/10k --slemm 10k --out 10k.chr$i.txt --chr $i; done

mv 10k.chr1.txt 10k.chrAll.txt
for i in `seq 2 30`; do tail -n +2 10k.chr$i.txt >> 10k.chrAll.txt; rm 10k.chr$i.txt; done
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
plink --r2 --bfile ../data/10k --chr-set 30 --out 10k --ld-window-r2 0 --ld-window 500 --ld-window-kb 10000
plink --assoc --bfile ../data/10k --allow-no-sex --chr-set 30 --pheno ../data/10k.bolt.txt --pheno-name QT --out 10k
```
