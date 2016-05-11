### Plotting and summarizing imputed metrics for only imputed variants either in
## 1) RefSeq gene regsion (i.e. exomic) OR
## 2) within 60kb of a custom SNPs

# this code generates a list of imputed variants in this target region
# this output will then be read into a custom version of ./Make_Final_Plots.R

#### SN 4/25/2014

rm(list=objects())
options(stringsAsFactors = FALSE, echo=TRUE)
dir <- "/projects/cidr/Ambrosone/"
library(rtracklayer)
library(IRanges)

## read in arguments - chrom and segment
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
seg <- args[2]

study <- "Ambrosone"

cat("Getting variant list for chrom",chr,"segment",seg,"\n")
   
## read in metrics for this segment
mets.fn <- paste(dir,"imputation/metrics/",study,"_chr",chr,".set",seg,".metrics",sep="")
mets <- read.table(mets.fn, header=TRUE)

## read in custom SNP list
cols <- rep("NULL", times=27)
cols[1:5] <- NA
cstm <- read.table(gzfile(paste(dir, "sarahcn/alleleMappings/20140411_check/Ambrosone_15041407_customExomeVariants_designManifestMerge.txt.gz",sep="")), header=T, colClasses=cols)

## read in segmentation scheme
segs <- read.csv(file=paste(dir,"imputation/dbGaP_doctn/imputation_segments.csv",sep=""))

## for table query, need "chr#" column
chr.char <- ifelse(segs$chrom!=23, segs$chrom, "X")
segs$chrom.ucsc <- paste("chr",chr.char, sep="")

## Use rtracklayer to access RefSeq records in this segment
session <- browserSession("UCSC")
genome(session) <- "hg19"
tmp <- segs[segs$chrom==chr & segs$segment==seg,]
grange <- GRanges(tmp$chrom.ucsc, IRanges(tmp$bp.start, tmp$bp.end))
cat("Querying RefSeq\n")
system.time(refSeq <- getTable(ucscTableQuery(session, track="RefSeq Genes", range=grange, table="refGene")))

# report query result
cat("\t",nrow(refSeq),"records found;",length(unique(refSeq$name2)),"unique gene names\n")

## get the position of type 2 custom SNP
cstm.type2 <- mets[is.element(mets$rs_id, cstm$Name) & mets$type==2,]
## add 60kb to either end
cstm.type2$start.pos <- cstm.type2$position-60000
cstm.type2$end.pos <- cstm.type2$position+60000
cat("Adding +/- 60 kb of",nrow(cstm.type2),"custom SNPs\n")
   
# define Target ranges
rangesTarget <- IRanges(start=c(refSeq$txStart, cstm.type2$start.pos),
                        end=c(refSeq$txEnd,cstm.type2$end.pos),
                        names=c(rep("RefSeq",times=nrow(refSeq)),
                          rep("CustomSNP",times=nrow(cstm.type2))))

# define Metrics ranges
rangesMets <- IRanges(mets$position, mets$position, names=mets$rs_id)

# query for overlap
ov <- overlapsAny(rangesMets, rangesTarget) 
table(ov)

# check that vector length is as expected
stopifnot(length(ov)==nrow(mets))

# flag type0 SNPs in overlap
mets$use <- mets$type==0 & ov

cat("Total of",sum(mets$use),"imputed variants in target regions, out of ",sum(mets$type==0),"total imputed\n")

## write out list of variants
var.write <- mets$rs_id[mets$use]

fn.out <- paste(dir,"sarahcn/plots/subsetPlot/variant_lists/vars_chr",chr,".seg",seg,".txt",sep="")
write.table(var.write, file=fn.out,col.names=F,row.names=F, quote=F,eol = "\n") 

