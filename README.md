# SLEMM
## Stochastic-Lanczos-Expedited Mixed Models

---

# Executable
https://github.com/jiang18/slemm/releases/tag/20220122

# Simulation
The procedures below generate genotypes and phenotypes of 500K individuals. Modify the sample size if needed. 
## Simulate genotypes with [genosim](https://aipl.arsusda.gov/software/genosim/)
1. pedigree.file and genotype.data0
```console
perl sim_ped.pl 500000
```
2. markersim.options
```
  50000      0      1.0         0           1
  loci     qtls    curvparm   correlated   traits

 96     30        .01       .99       0
seed  chromes  minfreq  maxfreq  QTLfreq
```
3. chips.txt
```
chip    reduce1   offset1    reduce2    offset2    depth1    error1
  1        1         0          1          0          35        .00
```
4. genosim.options
```
  1.0      .965      0
Morgans  LDparam  hapout

  30      30      .00      .000   .00000003
chromes Xchrome  notread   errate   mutate
```
5. Run genosim
```console
./markersim
./genosim
```
## Convert genosim output files to plink files
```console
perl aipl2plink.pl 500k
plink --file 500k --make-bed --out 500k --chr-set 30
```
## Simulate phenotypes based on genosim output files
```console
mkdir pheno
cd pheno
perl sim_snp_effects.pl 1.snp.csv
slemm --pred --binary_genotype ../500k --snp_estimate 1.snp.csv --output 1.gv
```
```console
Rscript --no-save sim_phe.R 1 0.3
```
# SLEMM
Note that the simulation procedures above generate a sample of 500K while the commands below work for a sample of 10K.
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
