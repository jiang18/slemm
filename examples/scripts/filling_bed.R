# This R script shows how to use bigsnpr to quickly fill missing genotypes in PLINK bed files.

# install.packages("bigsnpr")
library(bigsnpr)

bedfile <- "filename.bed"
filledbed <- "filename.filled.bed"

rds <- snp_readBed(bedfile, backingfile = tempfile())
bigsnp <- snp_attach(rds)
filled <- snp_fastImputeSimple(bigsnp$genotypes, method = "random", ncores = 10)
bigsnp$genotypes <- filled

snp_writeBed(bigsnp, filledbed)
