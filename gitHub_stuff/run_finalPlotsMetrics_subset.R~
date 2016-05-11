##########
## Plots and metrics
##########
# 4/18/2014

rm(list=objects())
options(stringsAsFactors = FALSE)
dir <- "/projects/cidr/Ambrosone/"

## preview some chrom's that have succesfully completed combining, for masked SNP test

imp_dir <- paste(dir,"imputation/", sep="")
out_dir <- paste(dir,"sarahcn/plots/", sep="")
study <- "Ambrosone"
schr <- 1
echr <- 23

source("/projects/geneva/geneva_sata/GCC_code/Imputation_IMPUTE2/R_functions/Make_Final_Plots.R")

## full set of graphs
make.finalimp.plot(imp_dir, out_dir, project=study, start_chr=schr, end_chr=echr)


## focus on masked SNPs to get strand check picture (i.e. sanity check on revised allele mappings)
# make.finalimp.plot(imp_dir, out_dir, project=study, start_chr=schr, end_chr=echr,
#                    sets="masked",plots="masked_check_strand")
