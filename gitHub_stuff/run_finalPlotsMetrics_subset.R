##########
## Plots and metrics
##########
# 4/28/2014

rm(list=objects())
options(stringsAsFactors = FALSE)
dir <- "/projects/cidr/Ambrosone/"

imp_dir <- paste(dir,"imputation/", sep="")

study <- "Ambrosone"
schr <- 1
echr <- 23

setwd <- "/projects/cidr/Ambrosone/sarahcn/R"

source("/projects/cidr/Ambrosone/sarahcn/R/Make_Final_Plots_Ambrosone.R")

#### Run once for imputed variants in exon targets
out_dir <- paste(dir,"sarahcn/plots/", sep="")
out_dir <- paste(dir,"sarahcn/plots/subsetPlot/exons/", sep="")
make.finalimp.plot(imp_dir, out_dir, project=study, start_chr=schr, end_chr=echr,
                   sets="imputed")

## #### Run once for imputed variants +/- 60kb of custom SNPs
## out_dir <- paste(dir,"sarahcn/plots/subsetPlot/customVars_60kbFlank/", sep="")
## make.finalimp.plot(imp_dir, out_dir, project=study, start_chr=schr, end_chr=echr,
##                    sets="imputed") 
