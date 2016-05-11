## From final imputation metrics reported by IMPUTE2, make standard plots for imputation report:
## (1) Quality metrics at all imputed SNPs
## (2) Quality metrics at all masked SNPs, from IMPUTE2 internal 'leave one out' masking expirement
## A combined version of Make_ImpMetrics_plots and Make_MaskedSNP_plots.R
## SN, 3/26/12

## 8/7/2012: Wenying Zheng revised masked SNP plot layout to have 3 plots
## 1/2/2013: SN revised to plot SNPs and indels/SVs separately for all imputed
## metrics (masked results still just report SNPs)
## 1/14/2014: SN streamlined code
## to loop through differnet variant sets (all, SNPs only, SVs and indels only)
## rather than repeat code
## 4/29/2014: SN added section to ancitipate metrics file
## dimensions based on qsub_combine_chr#.out file; also switching to using 1000G
## legend file to determine variant status (vs. grepping from .gprobs file, which
## was very slow and unwieldy)
## 6/30/2014: when reading in reference legend files, hard coding skipping over XTR region (in the future would be better to implement this by reading in imputation segmentation scheme)
## 8/24/2014: include indels in masked SNP summary (but not plots); also excluding from plots and summaries where r2_type0==-1 (i.e. undefined)
## 10/10/2014: for negative 'info' values, report number of "-1" observations and exclude them from the plots and summaries
## 12/12/2014: based on 1000G phase 3 legend files, changing alignment between metrics and refrence legend file to work in two parts: type 2 and then type 0 variants
## 12/30/2014: overhauled to split plots creation into several "sub functions"; also removing histograms
## 1/8/2015: changing legend file incorporation into per-chrom loop of reading in metrics file; vs reading all chroms metrics, then reading all chroms legend (to avoid alignment problems with SoL imputed 3)
## 1/20/2015: strand check graph - 1) report out # of masked SNPs in "problem quadrant" and 2) add density curves to both axes
## 2/11/2015: with IMPUTE2 v.2.3.2, a0 and a1 are in metrics files, so no longer need to load in IMPUTE2 legend files to determine variant alleles and type (i.e. SNP vs. indel)

### still to do: change from base graphics to ggplot2

#### Arguments:
## imp_dir -- directory name prefix in which 'metrics' folder can be found
## out_dir -- directory to write out plots and summaries
## sets -- defaults to processing both imputed and masked SNP results - otherwise can specify just one
## project -- name of project study dataset
## start_chr -- first chrom to include in plot (defaults to start at chr1)
## end_chr -- last chrom to include in plot (defaults to go through chrX)
## plots -- info_by_maf, and boxplot, to allow users to pick which plot types to create, for imputed and masked SNP metrics.
## fmt -- 'pdf' or 'png' format info by MAF; boxplots have to be png to enable %d option


#### Arguments specific to all imputed SNPs
## summary -- TRUE or FALSE - whether to print out SNP counts by chrom: all study, type 2, type 0+2

#### Arguments specific to masked SNPs
## info_filt -- value from [0,1], theoretical dosage r2 threshold for second column of graphs; default to 0.8
## verbose -- TRUE if you want to print out summary of mean and median masked SNP metrics for MAF<10% vs. MAF>=10%
## maf_thresh -- MAF value(s) used to dichotomize masked SNP results; default to 0.05

######## 


#### Define subfunctions

## I. Read in metrics files and annotate with variant type (SV/indel or SNP)
readMetrics <- function(imp_dir, project, start_chr, end_chr)

{
    message("Reading in imputation metrics on chroms ", start_chr, " through  ", 
        end_chr, "...")
    
    ## check that metrics files appear to be IMPUTE2 output
    hdr <- read.table(gzfile(paste0(imp_dir, "/metrics/combined/", project, "_chr", 
        start_chr, ".metrics.gz")), nrow = 1, header = TRUE, as.is = TRUE)
    missing_columns <- setdiff(c("rs_id", "info", "certainty", "type", "info_type0","concord_type0","r2_type0"), names(hdr))
    if (length(missing_columns) > 0) 
        stop("It does not look like you've given me IMPUTE2 metrics files - they're missing ", paste(missing_columns, collapse="; "))

    ## check for presence of a0 and a1 alleles
    missing_columns <- !is.element(c("a0", "a1"), names(hdr))
    if (sum(missing_columns) > 0) 
        stop("Starting with v2.3.2, IMPUTE2 metrics files contain a0 and a1. Your metrics files do not, so either re-run with IMPUTE2 v2.3.2 or later, or use an earlier (e.g., svn revision r1879) version of Make_Final_Plots.R -- which will utilize the reference legend file to determine variant alleles and type (SNP vs. indel).")
    
    ## use 'count.fields' to determine dimension of each per-chr metrics files and
    ## thus dimension of total combined chrom metrics file (rather than rbinding)
    dim.dat <- rep(NA, (end_chr - start_chr) + 1)
    names(dim.dat) <- start_chr:end_chr
    for (i in start_chr:end_chr) {
        mets.fn <- paste0(imp_dir, "/metrics/combined/", project, "_chr", i, ".metrics.gz")
        Cnt <- count.fields(mets.fn)
        nvars <- length(Cnt) - 1  # subtract header row
        dim.dat[names(dim.dat) == i] <- nvars
    }  # close chrom loop to determine variant dimensions

    ## loop through chosen chroms, combining metrics files
    ## add on column for chrom
    ## we'll add sv.indel and maf cols after looping through all chroms
    mets.all <- data.frame(matrix(nrow = sum(dim.dat), ncol = ncol(hdr) + 1))

    # start counter for keeping place in combined, annotated metrics files
    cnter <- 1
    for (i in start_chr:end_chr) {
        # report which chrom is being read in, in case of errors reading the combined
        # metrics file
        message("\tReading in metrics for chrom ", i, "...")
        nvars <- dim.dat[names(dim.dat) == i]
        mfil <- paste0(imp_dir, "/metrics/combined/", project, "_chr", i, ".metrics.gz")
        met <- read.table(gzfile(mfil), skip = 1, as.is = TRUE, comment.char = "", 
            nrow = nvars)
        # add chrom as final column
        met[, ncol(met) + 1] <- i
        names(met) <- c(names(hdr), "chr")

        ## add to other chroms' metrics files
        mets.all[cnter:(cnter + nvars - 1), ] <- met

        # use column names from first chrom
        if (i==start_chr) {colNamesFinal <- names(met)}

        # update counter
        cnter <- cnter + nvars
    } # close loop through chroms

    # use column names from first chrom
    names(mets.all) <- colNamesFinal

    ## Make variant type column: sv.indel or SNP
    mets.all$sv.indel <- FALSE
    # where either allele is > 1 nucleotide, assign to sv.indel=TRUE
    mets.all$sv.indel[nchar(mets.all$a0) > 1 | nchar(mets.all$a1) > 1] <- TRUE    

    ##### for genome-wide metrics file:
    message("\nTotal of ", sum(!mets.all$sv.indel), " SNPS and ", sum(mets.all$sv.indel), 
            " indels and structural variants...")

    ## check for negative info values - they can happen until Reed imputation in Oct
    ## 2014, hadn't observed them yet.  negative values indicate highly uncertain
    ## imputation; -1 means can't calculate info.  see
    ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#info_metric_details
    neg.info <- sum(mets.all$info < 0)
    neg1.info <- sum(mets.all$info == -1)
    if (neg.info > 0) {
        warning("There are ", neg.info, " variants with negative info score and ", 
            neg1.info, " variants with -1 (i.e. undefined) info scores. The undefined info value observations will not be included in the plots and summaries below, but should still be reported. Negative (but > -1) values will be included in the plots an summaries.\n")
    }
    
    ## Make MAF column
    mets.all$maf <- ifelse(mets.all$exp_freq_a1 > 0.5, 1 - mets.all$exp_freq_a1, 
        mets.all$exp_freq_a1)
    
    # check that there are at least 100 imputed variants
    nvar <- sum(mets.all$type == 0)
    if (nvar < 100) {
        stop("Come back when you have more imputed variants to graph, not just ", 
            nvar, "\n")
    }
    
    return(mets.all)
}


## II. Summarize imputed metrics
impSummary <- function(imp_dir, start_chr, end_chr, out_dir, project, mets.all) {
    message("Printing variant summary by chrom")
    ## total study snps read in initial SNP keep list, in case type 3 SNPs were
    ## excluded from IMPUTE2 output
    kept <- read.table(file = paste0(imp_dir, "/keeplists/snp.qualfilter.txt"), comment.char = "", 
        stringsAsFactors = FALSE)
    names(kept)[1:2] <- c("snp", "chr")
    ## count by chrom
    bychr <- as.data.frame(table(kept$chr))
    ## change chrom from factor to character
    bychr$Var1 <- as.character(bychr$Var1)
    ## convert chrom to numeric code
    bychr$Var1[bychr$Var1 == "X"] <- 23
    bychr$chrom.int <- as.numeric(bychr$Var1)
    bychr <- bychr[is.element(bychr$chrom.int, start_chr:end_chr), ]
    dat0 <- bychr[order(bychr$chrom.int), ]
    
    ## total input SNPs (imp basis + study only)
    dat1 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type, c(2, 3))]))
    ## total imputation basis SNPs
    dat2 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type, 2)]))
    ## total imputed SNPs
    dat3 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type, c(0, 2))]))
    
    dat1$study.snps.fmkeeplist <- dat0[, 2]
    dat1$imp.basis.snps <- dat2[, 2]
    dat1$imp.tot.snps <- dat3[, 2]
    names(dat1)[1:2] <- c("chrom", "study.snps")
    filo <- paste0(out_dir, project, "_chr", start_chr, "-", end_chr, ".SNPsummary.csv")
    write.csv(dat1, file = filo, eol = "\n", quote = FALSE, row.names = FALSE)
    
    # report fractions of imputed variants passing different info thresholds,
    # excluding info=-1 (i.e., undefined)
    met.imp <- mets.all[mets.all$type == 0 & mets.all$info > -1, ]
    message("\nSummaries of imputed variants passing various 'info' score thresholds:")
    
    ## loop through different info thresholds
    for (t in c(0.3, 0.5, 0.8)) {
        message("Fraction of SNPs with info >", t, " is ", round(sum(met.imp$info[!met.imp$sv.indel] > 
            t)/sum(!met.imp$sv.indel), 4))
        message("Fraction of SVs and indels with info >", t, " is ", round(sum(met.imp$info[met.imp$sv.indel] > 
            t)/sum(met.imp$sv.indel), 4))
        message("Fraction of all imputed variants with info >", t, " is ", round(sum(met.imp$info > 
            t)/nrow(met.imp), 4))
    }  # close loop on thresholds
}

## III. Graph imputed metrics: info by MAF
impPlotsInfo <- function(out_dir, start_chr, end_chr, fmt, maf, mets.ann) {
    # create datasets of SNPs and SVs/indels to plot separately
    mets.all.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1, ]
    mets.snp.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1 & !mets.ann$sv.indel, 
        ]
    mets.svi.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1 & mets.ann$sv.indel, 
        ]
    nsnps <- length(mets.snp.imp$rs_id)
    nsvi <- length(mets.svi.imp$rs_id)
    message("Info by MAF plots will represent a total of ", nsnps, " imputed SNPs and ", 
        nsvi, " imputed SVs and indels")
    
    # create output file name
    filo <- paste0(out_dir, "AllImpMetrics_chr", start_chr, "-", end_chr)

    mafs.save <- NULL
    
    ## loop through variant categories
    for (set in c("snp", "svs_indels", "all")) {
        if (set == "snp") 
            mets.curr <- mets.snp.imp
        if (set == "svs_indels") 
            mets.curr <- mets.svi.imp
        if (set == "all") 
            mets.curr <- mets.all.imp
        
        message("plotting info by MAF bin for ", set, " variants")

        ####### MAF data frame loop
        mafs.list <- seq(0, 0.5, length.out = 51)
        mafs <- as.data.frame(mafs.list[-51])
        names(mafs) <- "maf.min"
        mafs$nsnps <- mafs$mean.info <- rep(NA, times = length(mafs$maf.min))
        
        for (i in 1:length(mafs$maf.min)) {
            ## Count of variants per MAF bin
            mafs$nsnps[i] <- sum(mets.curr$maf > mafs.list[i] & mets.curr$maf <= 
                mafs.list[i + 1])
            ## Mean info
            mafs$mean.info[i] <- mean(mets.curr$info[mets.curr$maf > mafs.list[i] & 
                mets.curr$maf <= mafs.list[i + 1]])
        }
        ####### close MAF data frame loop
        
        ## plot values at the mid point of maf bin
        xval <- mafs$maf.min + 0.005
        
        if (fmt == "png") {
            png.fn <- paste(filo, set, "info_by_studyMAF.png", sep = ".")
            png(png.fn)
        }
        if (fmt == "pdf") {
            pdf.fn <- paste(filo, set, "info_by_studyMAF.pdf", sep = ".")
            pdf(pdf.fn)
        }
        
        ## Set margins. Default is $mar 5.1 4.1 4.1 2.1 -- increase the right margin for
        ## secondary y-axis
        par(mar = c(5.1, 4.1, 4.1, 4.1))
        plot(xval, mafs$mean.info, ylim = c(mafs$mean.info[1], 1), bty = "c", pch = 17, 
            xlab = "MAF of imputed variants", ylab = "info score (mean per MAF bin)")
        abline(h = seq(0, 1, by = 0.1), lty = 2, col = "darkgray")
        par(new = TRUE)
        plot(xval, mafs$nsnps, axes = FALSE, ylab = "", xlab = "", bty = "c", pch = 4, 
            col = "darkgreen")
        axis(4, tck = 0.01, las = 2)
        mtext("variant count", 4, col = "darkgreen")
        dev.off()
        
        # save data frame
        mafs$type <- set
        mafs.save <- rbind(mafs.save, mafs)
        
        
    }  # close loop on variant type
    
    # save MAF data frames to disk
    save(mafs.save, file = paste0(filo, "info_by_MAF.RData"))
}

## IV. Graph imputed metrics: by-chrom boxplots
impPlotsByChr <- function(out_dir, start_chr, end_chr, mets.ann) {
    # create datasets of SNPs and SVs/indels to plot separately
    mets.all.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1, ]
    mets.snp.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1 & !mets.ann$sv.indel, 
        ]
    mets.svi.imp <- mets.ann[mets.ann$type == 0 & mets.ann$info > -1 & mets.ann$sv.indel, 
        ]
    nsnps <- length(mets.snp.imp$rs_id)
    nsvi <- length(mets.svi.imp$rs_id)
    message("By-chrom box plots will represent a total of ", nsnps, " imputed SNPs and ", 
        nsvi, " imputed SVs and indels")
    
    # create output file name
    filo <- paste0(out_dir, "AllImpMetrics_chr", start_chr, "-", end_chr)
    
    ## loop through variant categories
    for (set in c("snp", "svs_indels", "all")) {
        
        if (set == "snp") 
            mets.curr <- mets.snp.imp
        if (set == "svs_indels") 
            mets.curr <- mets.svi.imp
        if (set == "all") 
            mets.curr <- mets.all.imp
        
        message("making boxplots for ", set, " variants")
        png.fn <- paste(filo, set, "boxplot_%d.png", sep = ".")
        png(png.fn)
        
        boxplot(info ~ chr, data = mets.curr, outline = FALSE, xlab = "Chromosome", 
            ylab = "info", main = NULL, las = 3, col = "seashell2", cex.axis = 0.9)
        
        boxplot(certainty ~ chr, data = mets.curr, outline = FALSE, xlab = "Chromosome", 
            ylab = "certainty", main = NULL, las = 3, col = "seashell2", cex.axis = 0.9)
        
        boxplot(maf ~ chr, data = mets.curr, outline = FALSE, xlab = "Chromosome", 
            ylab = "imputed MAF", main = NULL, las = 3, col = "seashell2", cex.axis = 0.9)
        
        dev.off()
        
    }  # close loop through variant types
}

## V. Graph masked metrics
maskedPlots <- function(out_dir, start_chr, end_chr, maf, info_filt, mets.ann) {
    ## Count # of SNPs in internal masked SNP expirement
    mets.mask <- mets.ann[mets.ann$type == 2 & !mets.ann$sv.indel, ]
    nsnps <- length(mets.mask$rs_id)
    message("\n\nTotal number of masked SNPs will be ", nsnps)
    
    # report number of masked SVs or indels
    mets.svi.mask <- mets.ann[mets.ann$type == 2 & mets.ann$sv.indel, ]
    nsvi <- length(mets.svi.mask$rs_id)
    message("\tThere are a total of ", nsvi, " masked SVs and indels, which will be excluded from plots but included in the summary table...")
    
    # create base name for all masked metrics plots
    filo <- paste(out_dir, "AllMaskedSNPMetrics_chr", start_chr, "-", end_chr, sep = "")
    use <- mets.mask
    ####### MAF data frame loop use function to loop through each MAF bin and find mean
    ####### masked SNP metrics values
    mafs.list <- seq(0, 0.5, length.out = 51)
    mafs <- as.data.frame(mafs.list[-51])
    names(mafs) <- "maf.min"
    mafs$nsnps <- mafs$meanr2 <- mafs$meanc <- rep(NA, times = length(mafs$maf.min))
    
    for (i in 1:length(mafs$maf.min)) {
        ## Count of SNPs per MAF bin
        mafs$nsnps[i] <- sum(use$maf > mafs.list[i] & use$maf <= mafs.list[i + 1])
        ## Mean dosage r2 - exclude undefined r2_type0 (i.e. r2_type0=-1)
        mafs$meanr2[i] <- mean(use$r2_type0[use$r2_type0 != -1 & use$maf > mafs.list[i] & 
            use$maf <= mafs.list[i + 1]])
        ## Mean concordance
        mafs$meanc[i] <- mean(use$concord_type0[use$maf > mafs.list[i] & use$maf <= 
            mafs.list[i + 1]])
    }
    ####### MAF data frame loop
    
    ## remove MAF bins with no snps
    mafs <- mafs[mafs$nsnps > 0, ]
    mafs.all <- mafs
    
    ## repeat thresholding on info_type0
    use.r2 <- mets.mask[mets.mask$info_type0 >= info_filt, ]
    dim(use.r2)
    
    ####### MAF data frame loop use function to loop through each MAF bin and find mean
    ####### concordance values
    mafs.list <- seq(0, 0.5, length.out = 51)
    mafs <- as.data.frame(mafs.list[-51])
    names(mafs) <- "maf.min"
    mafs$nsnps <- mafs$meanr2 <- mafs$meanc <- rep(NA, times = length(mafs$maf.min))
    
    for (i in 1:length(mafs$maf.min)) {
        ## Count of SNPs per MAF bin
        mafs$nsnps[i] <- sum(use.r2$maf > mafs.list[i] & use.r2$maf <= mafs.list[i + 
            1])
        ## Mean dosage r2 - exclude undefined r2_type0 (i.e. r2_type0=-1)
        mafs$meanr2[i] <- mean(use.r2$r2_type0[use.r2$r2_type0 != -1 & use.r2$maf > 
            mafs.list[i] & use.r2$maf <= mafs.list[i + 1]])
        ## Mean concordance
        mafs$meanc[i] <- mean(use.r2$concord_type0[use.r2$maf > mafs.list[i] & use.r2$maf <= 
            mafs.list[i + 1]])
    }
    ####### MAF data frame loop
    
    mafs.dosr2 <- mafs
    mafs.dosr2$pct.snps <- round(mafs.dosr2$nsnps/mafs.all$nsnps, 3)
    
    ## combine two data frames of average metrics by MAF bin
    mafs.all$pct.snps <- 1
    mafs.comb <- rbind(mafs.all, mafs.dosr2)
    save(mafs.comb, file = paste(filo, "mets_by_MAF.RData", sep = ""))
    
    ## plot values at the mid point of maf bin
    xval <- mafs$maf.min + 0.005
    
    #### masked SNP plots: (1) SNP count; (2) Dosage; (3) Concordance with and without
    #### an info_type0 threshold 1x3 landscape plot (per Wenying Zheng revision,
    #### 8/7/2012)
    
    pdf.fn <- paste(filo, "comboplot.pdf", sep = ".")
    pdf(pdf.fn, width = 20, height = 8, paper = "a4r")
    par(mfcol = c(1, 3), mar = c(1, 4, 1, 3), pty = "s")  ## version for square
    
    ## plot for snp count
    plot(x = xval, y = mafs.all$nsnps, bty = "c", pch = 16, xlab = "", ylab = "", 
        lwd = 2)
    mtext("SNP Count", 2, line = 2, cex = 0.8)
    mtext("Study MAF", 1, line = 2.5, cex = 0.8, font = 3)
    par(new = TRUE)
    plot(xval, mafs.dosr2$pct.snps * 100, axes = FALSE, ylab = "", xlab = "", bty = "c", 
        pch = 16, col = "darkgray", lwd = 2)
    axis(4, at = c(0, 110), labels = c("", ""), lwd.ticks = 0)
    axis(4)
    mtext("% SNP from bin", 4, line = 2, cex = 0.8, col = "darkgray")
    legend("right", c("all SNPs", paste("SNPs with info_type0>=", info_filt)), pch = c(16, 
        16), col = c("black", "darkgray"), bg = "white")
    
    ## plot for dosage
    plot(x = xval, y = mafs.all$meanr2, ylim = c(min(mafs.all$meanr2, na.rm = TRUE), 
        1), pch = 16, xlab = "", ylab = "")
    mtext("Dosage r2", 2, line = 2, cex = 0.8)
    mtext("Study MAF", 1, line = 2, cex = 0.8, font = 3)
    abline(h = seq(0, 1, by = 0.1), lty = 2, col = "darkgray")
    points(x = xval, y = mafs.dosr2$meanr2, pch = 16, col = "darkgray")
    legend("bottomright", c("all SNPs", paste("SNPs with info_type0>=", info_filt)), 
        pch = c(16, 16), col = c("black", "darkgray"), bg = "white")
    
    ## plot for concordance
    plot(x = xval, y = mafs.all$meanc, ylim = c(min(mafs.all$meanc, na.rm = TRUE), 
        1), pch = 16, xlab = "", ylab = "")
    mtext("Concordance", 2, line = 2, cex = 0.8)
    mtext("Study MAF", 1, line = 2, cex = 0.8, font = 3)
    abline(h = seq(0, 1, by = 0.01), lty = 2, col = "darkgray")
    points(x = xval, y = mafs.dosr2$meanc, pch = 16, col = "darkgray")
    legend("topright", c("all SNPs", paste("SNPs with info_type0>=", info_filt)), 
        pch = c(16, 16), col = c("black", "darkgray"), bg = "white")
    
    dev.off()
}

## VI. Summarize masked metrics
maskedSummary <- function(maf_thresh, mets.ann) {
    mets.snp.mask <- mets.ann[mets.ann$type == 2 & !mets.ann$sv.indel, ]
    mets.svi.mask <- mets.ann[mets.ann$type == 2 & mets.ann$sv.indel, ]
    
    for (m in maf_thresh) {
        # loop through SNPs and then SV/indels
        for (type in c("SNP", "SV/INDEL")) {
            message("Dichotomizing masked ", type, " metrics by study MAF of ", maf_thresh)
            
            if (type == "SNP") {
                mets.mask <- mets.snp.mask
            }
            if (type == "SV/INDEL") {
                mets.mask <- mets.svi.mask
            }
            
            ## Report mean and median concordance
            message("Summary of concord_type0")
            mna <- round(mean(mets.mask$concord_type0[mets.mask$maf < maf_thresh]), 
                4)
            mda <- round(median(mets.mask$concord_type0[mets.mask$maf < maf_thresh]), 
                4)
            mnb <- round(mean(mets.mask$concord_type0[mets.mask$maf >= maf_thresh]), 
                4)
            mdb <- round(median(mets.mask$concord_type0[mets.mask$maf >= maf_thresh]), 
                4)
            
            message(">>>>> MAF < ", maf_thresh, ", mean=", mna, " median=", mda, 
                "...\n")
            message(">>>>> MAF >=", maf_thresh, ", mean=", mnb, " median=", mdb, 
                "...\n")
            
            ## Report mean and median dosage r2 (r2_type0)
            ## Exclude where r2_type0 is undefined
            nmono.missR2 <- sum(mets.mask$r2_type0 == -1 & mets.mask$maf == 0)
            ## undefined r2_type0 should only happen at monomorphic variants,
            ## ...but confirm and issue a warning if non-monomorphic variants have undefined r2_type0
            nonmono.missR2 <- sum(mets.mask$r2_type0 == -1 & mets.mask$maf > 0)
            if (nonmono.missR2 > 0) {
                warning("There are ", nonmono.missR2, " masked ", type, " variants that have r2_type0=-1, despite being polymorphic. Note these are excluded from the plots and summaries of r2_type0")
            }
            message("Summary of r2_type0, excluding ", nmono.missR2, " monomorphic variants with undefined r2_type0")
            
            mets.use <- mets.mask[mets.mask$r2_type0 != -1, ]
            mna <- round(mean(mets.use$r2_type0[mets.use$maf < maf_thresh]), 4)
            mda <- round(median(mets.use$r2_type0[mets.use$maf < maf_thresh]), 4)
            mnb <- round(mean(mets.use$r2_type0[mets.use$maf >= maf_thresh]), 4)
            mdb <- round(median(mets.use$r2_type0[mets.use$maf >= maf_thresh]), 4)
            
            message(">>>>> MAF < ", maf_thresh, ", mean=", mna, " median=", mda, 
                "...\n")
            message(">>>>> MAF >= ", maf_thresh, ", mean=", mnb, " median=", mdb, 
                "...\n")
            
            nlowmaf <- sum(mets.mask$maf < maf_thresh)
            nhimaf <- sum(mets.mask$maf >= maf_thresh)
            nmono <- sum(mets.mask$maf == 0)
            message("Count w/study MAF < ", maf_thresh, ": ", nlowmaf, "\n")
            message("Count w/study MAF >= ", maf_thresh, ": ", nhimaf, "\n")
            message("Count of monomorphic: ", nmono, "\n\n")
            
        }  ## close loop on variant type (SNP vs. SV/indel)
    }  ## close loop through MAF thresholds
}

## VII. Masked strand check
maskedStrandCheck <- function(out_dir, start_chr, end_chr, mets.ann) {
    filo <- paste(out_dir, "AllMaskedSNPMetrics_chr", start_chr, "-", end_chr, sep = "")
    
    png.fn <- paste(filo, "check_strand.png", sep = ".")
    
    # limit to masked SNPs (exclude masked SVs/indels) with non-missing info_type0
    mets.mask <- mets.ann[mets.ann$type == 2 & !mets.ann$sv.indel &
                          mets.ann$info_type0 >= 0, ]
    
    # flag strand ambiguous SNPs
    mets.mask$alleles <- pasteSorted(mets.mask$a0, mets.mask$a1)
    mets.mask$strand.amb <- is.element(mets.mask$alleles, c("A/T", "C/G"))

    ## ggplot: High Density Scatterplot with Color Transparency
    cols.transp <- rgb(0, 100, 0, 50, maxColorValue = 255)
    leg.name <- "Strand Ambiguous\n(A/T or C/G) SNPs"
    # make scatterplot
    plot_scatter <- ggplot(mets.mask, aes(x = info_type0, y = concord_type0,
                          color=strand.amb, shape=strand.amb)) +
                    geom_point() + scale_shape_manual(values=c(16,8), name=leg.name) +
                    scale_color_manual(values=c(cols.transp,"black"), name=leg.name) +
                    xlim(0,1) + ylim(0,1) + theme(legend.position="bottom") +
                    geom_abline(intercept = 0, slope = 1, linetype="dashed")

    # see http://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/
    #placeholder plot - prints nothing at all
    empty <- ggplot()+geom_point(aes(1,1), colour="white") +
         theme(                              
           plot.background = element_blank(), 
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
           panel.border = element_blank(), 
           panel.background = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank()
         )    

    #marginal density of x - plot on top
    plot_top <- ggplot(mets.mask, aes(info_type0)) +
      geom_density() + theme(legend.position = "none") +
        xlab("") + ylab("")

    #marginal density of y - plot on the right
    plot_right <- ggplot(mets.mask, aes(concord_type0)) + geom_density() +
      coord_flip() + theme(legend.position = "none") +
        xlab("") + ylab("")

    #arrange the plots together, with appropriate height and width for each row and column
    # note ggsave doesn't work with grid.arrange
    png(png.fn)
    grid.arrange(plot_top, empty, plot_scatter, plot_right, ncol=2, nrow=2,
                 widths=c(4, 1), heights=c(1, 4))
    dev.off()
  
    # report out number of type 2 SNPs in "problem quadrant"
    nprob <- sum(mets.mask$concord_type0<0.5 & mets.mask$info_type0>0.5)
    nprob.amb <- sum(mets.mask$concord_type0<0.5 &
                     mets.mask$info_type0>0.5 & mets.mask$strand.amb)
    message("Of ",nrow(mets.mask)," total masked SNPs, ",nprob," are in problem quadrant of strand check (info_type0>0.5 and concord_type0<0.5) - ", nprob.amb, " of which are also strand ambiguous SNPs")
}

######## End subfunctions
######## Define main function (wrapper function that calls subfunctions)

make.finalimp.plot <- function(imp_dir, out_dir, project, start_chr = 1, end_chr = 23, 
    sets = c("imputed", "masked"), plots = c("boxplot", "info_by_maf", "masked_check_strand"), 
    fmt = "pdf", summary = TRUE, info_filt = 0.8, verbose = TRUE, maf_thresh = 0.05) 
{
    
    # should not have to explicitly load these libraries, but in the past have gotten
    # sporadic 'cannot find function' messages
    require("graphics")
    require("stats")
    # require('gdsfmt', lib.loc = '/projects/geneva/gcc-fs2/R_packages/library/')
    require("GWASTools")  #, lib.loc = '/projects/geneva/gcc-fs2/R_packages/library/')
    require("ggplot2")
    require("gridExtra")
    options(stringsAsFactors = FALSE)
    
    # I. Read in metrics and annotate with alleles
    mets.annotated <- readMetrics(imp_dir, project, start_chr, end_chr)
    
    ## II. Summarize imputed metrics
    if (is.element("imputed", sets) & summary) {
        impSummary(imp_dir, start_chr, end_chr, out_dir, project, mets.all = mets.annotated)
    }
    
    ## III. Graph imputed metrics: info by MAF
    if (is.element("imputed", sets) & is.element("info_by_maf", plots)) {
        impPlotsInfo(out_dir, start_chr, end_chr, fmt, maf = maf_thresh, mets.ann = mets.annotated)
    }
    
    ## IV. Graph imputed metrics: by-chrom boxplots
    if (is.element("imputed", sets) & is.element("boxplot", plots)) {
        impPlotsByChr(out_dir, start_chr, end_chr, mets.ann = mets.annotated)
    }
    
    ## V. Graph masked metrics
    if (is.element("masked", sets)) {
        maskedPlots(out_dir, start_chr, end_chr, maf = maf_thresh, info_filt, mets.ann = mets.annotated)
    }
    
    ## VI. Summarize masked metrics
    if (is.element("masked", sets) & verbose) {
        maskedSummary(maf_thresh, mets.ann = mets.annotated)
    }
    
    ## VII. Masked strand check
    if (is.element("masked", sets) & is.element("masked_check_strand", plots)) {
        maskedStrandCheck(out_dir, start_chr, end_chr, mets.ann = mets.annotated)
    }
    
}  # close main wrapper function     
