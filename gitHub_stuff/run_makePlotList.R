### Plotting and summarizing imputed metrics for only imputed variants
## Making two different lists
## 1) Exons (from RefSeq gene records)
## 2) within 60kb of a custom SNPs

# this code generates a list of imputed variants in these target regions

## [sarahcn@pearson0 subsetPlot]$ pwd
## /projects/cidr/Ambrosone/sarahcn/plots/subsetPlot
## [sarahcn@pearson0 subsetPlot]$ lsd
## drwxr-sr-x+ 3 sarahcn cidr 3 Apr 28 09:33 customVars_60kbFlank
## drwxr-sr-x+ 3 sarahcn cidr 3 Apr 28 09:33 exons

# this output will then be read into a custom version of ./Make_Final_Plots.R

#### SN 4/28/2014

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

# define Metrics ranges
# ...only for type0 SNPs
mets.imp <- mets[mets$type==0,]
rangesMets <- IRanges(mets.imp$position, mets.imp$position, names=mets.imp$rs_id)

## read in segmentation scheme
segs <- read.csv(file=paste(dir,"imputation/dbGaP_doctn/imputation_segments.csv",sep=""))

## for table query, need "chr#" column
chr.char <- ifelse(segs$chrom!=23, segs$chrom, "X")
segs$chrom.ucsc <- paste("chr",chr.char, sep="")


################# I. Query and write out imputed variants in exons

## Use rtracklayer to access RefSeq records in this segment
session <- browserSession("UCSC")
genome(session) <- "hg19"
tmp <- segs[segs$chrom==chr & segs$segment==seg,]
grange <- GRanges(tmp$chrom.ucsc, IRanges(tmp$bp.start, tmp$bp.end))
cat("Querying RefSeq\n")
system.time(refSeq <- getTable(ucscTableQuery(session, track="RefSeq Genes", range=grange, table="refGene")))

# report query result
cat("\t",nrow(refSeq),"records found;",length(unique(refSeq$name2)),"unique gene names\n")

# extract exon start positions and exon end positions
# head(refSeq$exonStarts)

# extract all the exon start coordinates into one column
list.start <- strsplit(refSeq$exonStarts, ",")
matx.starts <- matrix(unlist(list.start), ncol=1, byrow=TRUE)

# extract all the exon end coordinates into one column
list.end <- strsplit(refSeq$exonEnds, ",")
matx.ends <- matrix(unlist(list.end), ncol=1, byrow=TRUE)

# combine columns and make numeric
matx.exons <- cbind(as.numeric(matx.starts), as.numeric(matx.ends))

# check that dimensions of this matx exons matches expected number of exons
exon.cnt.expected <- sum(refSeq$exonCount)

if(nrow(matx.exons)!=exon.cnt.expected)
     {stop("Number of extracted exons coordinates (",nrow(matx.exons),
           ") does not match sum of exonCount column (", exon.cnt.expected,")")}

# define Target ranges (all exons)
rangesTarget <- IRanges(start=matx.exons[,1], end=matx.exons[,2])

# query for overlap with imptued Variants
ov <- overlapsAny(rangesMets, rangesTarget) 
table(ov)

# check that vector length is as expected
stopifnot(length(ov)==nrow(mets.imp))

cat("Total of",sum(ov),"imputed variants in exon target regions, out of",sum(mets$type==0),"total imputed\n")

## write out list of variants to use: type0 SNPs in overlap
var.write <- mets.imp$rs_id[ov]

fn.out <- paste(dir,"sarahcn/plots/subsetPlot/exons/variant_lists/vars_chr",chr,".seg",seg,".txt",sep="")
write.table(var.write, file=fn.out,col.names=F,row.names=F, quote=F,eol = "\n") 

################# II. Query and write out imputed variants with 60kb of custom variants

## read in custom SNP list
cols <- rep("NULL", times=27)
cols[1:5] <- NA
cstm <- read.table(gzfile(paste(dir, "sarahcn/alleleMappings/20140411_check/Ambrosone_15041407_customExomeVariants_designManifestMerge.txt.gz",sep="")), header=T, colClasses=cols)

## get the position of type 2 custom SNP
cstm.type2 <- mets[is.element(mets$rs_id, cstm$Name) & mets$type==2,]

## add 60kb to either end
cstm.type2$start.pos <- cstm.type2$position-60000
cstm.type2$end.pos <- cstm.type2$position+60000
cat("Adding +/- 60 kb of",nrow(cstm.type2),"custom SNPs\n")
   
# define Target ranges (+/- 60kb of custom SNPs)
# per Cathy: "There should still be some LD at about 60 kb; in most regions of the genome it occurs within about 100 kb"
rangesTarget <- IRanges(start=cstm.type2$start.pos,end=cstm.type2$end.pos)

# query for overlap
ov <- overlapsAny(rangesMets, rangesTarget) 
table(ov)

# check that vector length is as expected
stopifnot(length(ov)==nrow(mets.imp))

cat("Total of",sum(ov),"imputed variants w/in 60kb of a custom SNP, out of",sum(mets$type==0),"total imputed\n")

## write out list of variants to use: type 0 SNPs within 60kb of custom SNP
var.write <- mets.imp$rs_id[ov]

fn.out <- paste(dir,"sarahcn/plots/subsetPlot/customVars_60kbFlank/variant_lists/vars_chr",chr,".seg",seg,".txt",sep="")
write.table(var.write, file=fn.out,col.names=F,row.names=F, quote=F,eol = "\n") 
