# Simulation
## Simulate genotypes with [genosim](https://aipl.arsusda.gov/software/genosim/)
```console
cd data
perl ../scripts/sim_ped.pl 10000
markersim
genosim
```
- Enter the data folder.
- The Perl script generates pedigree.file and genotype.data0 for 10k unrelated individuals. 
- Run genosim.
## Convert genosim output files to PLINK files
```console
perl ../scripts/aipl2plink.pl 10k
plink --file 10k --make-bed --out 10k --chr-set 30
```
## Simulate phenotypes based on genotypes
```console
perl ../scripts/sim_snp_effects.pl 1.snp.csv
slemm --pred --binary_genotype 10k --snp_estimate 1.snp.csv --output 10k.gv.csv
Rscript --no-save ../scripts/sim_phe.R 10k 0.3
```
- The Perl script simulates SNP effects.
- slemm computes total genetic values.
- The R script simulates phenotypes with a heritability of 0.3.
# SLEMM
## SNP info file
```console
cd ..
mkdir slemm
cd slemm
perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < ../data/10k.bim > snp_info.csv
```
- The SNP info file is a CSV file with a header line and specifies what SNPs to be included in SLEMM.
## REML
```console
slemm --reml --phenotype_file ../data/10k.slemm.csv --binary_genotype_file ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
```
## Prediction of genomic breeding values
```console
slemm --pred --binary_genotype_file ../data/10k --snp_estimate 10k.reml.snp.csv --out 10k.gebv.csv
```
## LMM and GWA
```console
slemm --lmm --phenotype_file ../data/10k.slemm.csv --binary_genotype_file ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
OMP_NUM_THREADS=1 slemm_gamma.py --pfile ../data/10k --slemm 10k --out 10.gamma.txt
OMP_NUM_THREADS=10 slemm_gwa.py --pfile ../data/10k --slemm 10k --out 10k.chr1.txt --chr 1
```
- slemm --lmm fits linear mixed model for genome-wide associations.
- slemm_gamma.py computes GRAMMAR-Gamma association statistics for each SNP using the output of --lmm.
- slemm_gwa.py computes single-SNP association statistics closely approximating EMMAX. 

```bash
export OMP_NUM_THREADS=10
for i in `seq 1 30`; do slemm_gwa.py --pfile ../data/10k --slemm 10k --out 10k.chr$i.txt --chr $i; done

cp 10k.chr1.txt 10k.chrAll.txt
for i in `seq 2 30`; do tail -n +2 10k.chr$i.txt >> 10k.chrAll.txt; done
```
- Association tests for all chromosomes
# [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
```console
mkdir bolt
cd bolt

~/software/BOLT-LMM_v2.3.4/bolt \
  --lmmInfOnly --bed=../10k/10k.bed --bim=../10k/10k.bim --fam=../10k/10k.fam \
  --phenoFile=../10k/pheno/1.bolt.txt --phenoCol=QT \
  --numThreads=10 --Nautosomes=30 --LDscoresUseChip \
  --verboseStats --statsFile=10k.1.lmm.txt
```
# [PLINK](https://www.cog-genomics.org/plink/1.9/)
```console
mkdir plink
cd plink
plink --r2 --bfile ../10k/10k --chr-set 30 --out 1 --ld-window-r2 0 --ld-window 500 --ld-window-kb 10000
plink --assoc --bfile ../10k/10k --allow-no-sex --chr-set 30 --pheno ../10k/pheno/1.bolt.txt --pheno-name QT --out 1
```
