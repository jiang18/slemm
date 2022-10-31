# Simulation
## Simulate genotypes with [genosim](https://aipl.arsusda.gov/software/genosim/)
1. Generate pedigree.file and genotype.data0 for genosim
```console
cd data
perl ../scripts/sim_ped.pl 10000
```
The procedure above generates genosim input files for 10k unrelated individuals. 
2. Run genosim
```console
markersim
genosim
```
## Convert genosim output files to plink files
```console
perl ../scripts/aipl2plink.pl 10k
plink --file 10k --make-bed --out 10k --chr-set 30
```
## Simulate phenotypes based on genosim output files
```console
perl ../scripts/sim_snp_effects.pl 1.snp.csv
slemm --pred --binary_genotype 10k --snp_estimate 1.snp.csv --output 10k.1.gv.csv
Rscript --no-save sim_phe.R 10k.1 0.3
```
- The Perl script simulates SNP effects.
- slemm computes total genetic values.
- The R script simulates phenotypes with a heritability of 0.3.
# SLEMM
## REML
```console
mkdir slemm
cd slemm
slemm --phenotype_file ../10k/pheno/1.slemm.csv --binary_genotype_file ../10k/10k --trait QT --reml --snp_info_file snp_info.csv --out 1 --num_threads 10 --num_random_probes 30
```
## Weighted least squares
```console
slemm --phenotype_file ../10k/pheno/1.slemm.csv --binary_genotype_file ../10k/10k --trait QT --wls --snp_info_file snp_info.csv --out 1 --num_threads 10
```
## LMM and GWA
```console
slemm --phenotype_file ../10k/pheno/1.slemm.csv --binary_genotype_file ../10k/10k --trait QT --lmm --snp_info_file snp_info.csv --out 1 --num_threads 10 --num_qf_markers 30
slemm_gamma.py --pfile ../10k/10k --slemm 1 --out 1.gamma.txt
OMP_NUM_THREADS=10 slemm_gwa.py --pfile ../10k/10k --slemm 1 --out 1.chr1.txt --chr 1
```

Association tests for chromosomes 1-22
```bash
export OMP_NUM_THREADS=10
for i in `seq 1 22`; do slemm_gwa.py --pfile ../10k/10k --slemm 1 --out 1.chr$i.txt --chr $i; done

cp 1.chr1.txt 1.chr1-22.txt
for i in `seq 2 22`; do tail -n +2 1.chr$i.txt >> 1.chr1-22.txt; done
```
## Dependencies of slemm_gamma.py and slemm_gwa.py
- Linux packages: python3, python3-devel, gcc, and gcc-c++
- Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)

It is usually easy to fullfill the dependencies if your machine has python3 and pip. If it is difficult for you to pip the python3 packages, please try the ones I have compiled on CentOS Linux 7 with GCC 4.8.5 (https://github.com/jiang18/slemm#executable). 

## Error weight
The argument ```--error_weight_name``` can be added in --reml/--wls/--lmm to weight individual error variance. To enable this, add one additional column in the phenotype CSV file and specify the header name by --error_weight_name. Values are typically 1/r^2-1, where r^2 is the reliability of deregressed EBVs. 
# [GCTA-fastGWA](https://cnsgenomics.com/software/gcta/#fastGWA)
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
