############ Custom version for Ambrosone to limit to exome and other targeted variants

## From final imputation metrics reported by IMPUTE2, make standard plots for imputation report:
	## (1) Quality metrics at all imputed SNPs
	## (2) Quality metrics at all masked SNPs, from IMPUTE2 internal 'leave one out' masking expirement
## A combined version of Make_ImpMetrics_plots and Make_MaskedSNP_plots.R
## SN, 3/26/12
## 8/7/2012: Wenying Zheng revised masked SNP plot layout to have 3 plots
## 1/2/2013: SN revised to plot SNPs and indels/SVs separately for all imputed metrics (masked results still just report SNPs)
## 1/14/2014: SN streamlined code too loop through differnet variant sets (all, SNPs only, SVs and indels only) rather than repeat code

## Arguments:
## imp_dir -- directory name prefix in which "metrics" folder can be found
## out_dir -- directory to write out plots and summaries
## sets -- defaults to processing both imputed and masked SNP results - otherwise can specify just one
## project -- name of project study dataset
## start_chr -- first chrom to include in plot (defaults to start at chr1)
## end_chr -- last chrom to include in plot (defaults to go through chrX)
## plots -- hist, info_by_maf, and boxplot, to allow users to pick which of three plot types, for imputed and masked SNP metrics.
        ## as of Jan 2013, removing "hist" from default list of plot types, as the histograms of info and certainty across all SNPs aren't usually that informative
        ## "masked_check_strand" will plot concordance vs info (impuation confidence) for masked SNPs, to assess possible strand issues
## fmt -- "pdf" or "png" format for histogram and info by MAF; boxplots have to be png to enable %d option

#### Arguments specific to all imputed SNPs SNPs
## summary -- TRUE or FALSE - whether to print out SNP counts by chrom: all study, type 2, type 0+2

#### Arguments specific to masked SNPs
## info_filt -- value from [0,1], theoretical dosage r2 threshold for second column of graphs; default to 0.8
## verbose -- TRUE if you want to print out summary of mean and median masked SNP metrics for MAF<10% vs. MAF>=10%
## maf_thresh - MAF value(s) used to dichotomize masked SNP results; default to 0.05

########

make.finalimp.plot <- function(imp_dir, out_dir, project, start_chr=1, end_chr=23, sets=c("imputed","masked"),
							plots=c("boxplot","info_by_maf", "masked_check_strand"), fmt="pdf", summary=TRUE,
							info_filt=0.8, verbose=TRUE, maf_thresh=0.05)
  {

  # should not have to explicitly load these libraries, but in the past have gotten sporadic "cannot find function" messages
    require("graphics")
    require("stats")
    require("ggplot2")
    options(stringsAsFactors = FALSE)
    
    cat("Reading in imputation metrics on chroms", start_chr, "through", end_chr, "...\n")
    
    ## check that metrics files appear to be IMPUTE2 output
    hdr <- read.table(gzfile(paste(imp_dir,"/metrics/combined/",project ,"_chr",start_chr,".metrics.gz",sep="")),nrow=1, header=TRUE, as.is=TRUE)
    missing_columns <- !is.element(c("rs_id", "info", "certainty", "type"), names(hdr))
    if(sum(missing_columns) > 0)
       stop("It does not look like you've given me IMPUTE2 metrics files")

    ## determine dimension of each per-chr metrics files and thus dimension of total combined chrom metrics file (rather than rbinding)
    # parse through log file of combining metrics to get variant count
    dim.dat <- rep(NA, (end_chr-start_chr)+1)
    names(dim.dat) <- start_chr:end_chr
    for (i in start_chr:end_chr){
    cmnd <- paste("awk '{print NR, $0}'", paste(imp_dir,"/imputed/combined/qsub_combine_",i,".out",sep=""),
                  "| egrep 'combined metrics' | awk '{print $1}'")
    getLine <- as.integer(system(cmnd,intern=TRUE))

    # the next line contains the variant cound
    cmnd <- paste("awk 'NR==",getLine+1," {print $0}' ", paste(imp_dir,"/imputed/combined/qsub_combine_",i,".out",sep=""),sep="")
    getCount <- as.integer(system(cmnd,intern=TRUE))

    dim.dat[names(dim.dat)==i] <- getCount

  } # close chrom loop to determine variant dimensions                                         

    ## loop through chosen chroms, combining metrics files
    mets.all <- data.frame(matrix(nrow=sum(dim.dat), ncol=ncol(hdr)+1))

    # start counter for keeping place in combined metrics files
    cnter <- 1
    
    for (i in start_chr:end_chr){

      # report which chrom is being read in, in case of errors reading the combined metrics file
      cat("Reading chrom", i, "...\n")

      nvars <- dim.dat[names(dim.dat)==i]
      
      mfil <- paste(imp_dir,"/metrics/combined/",project ,"_chr",i,".metrics.gz",sep="")
      met <- read.table(gzfile(mfil), skip=1, as.is=TRUE, comment.char="",
                        nrow=nvars)
      # add chrom as final column
      met[,ncol(met)+1] <- i
      mets.all[cnter:(cnter+nvars-1),] <- met

      # update counter
      cnter <- cnter+nvars
    }

    names(mets.all) <- c(names(hdr), "chr")

    ## check for negative info values - they can happen, but I've not observed them yet and thus this plotting and summarizing may not deal with them appropriately. If there are negative info values, issue a warning and continue plotting
    neg.info <- sum(mets.all$info<0)
    if(neg.info>0) warning("There are ",neg.info," variants with negative info score; these plots and summaries may not deal with these appropriately.\n I suggest revisiting and possibly modifying this function.\n")

    ########## for Ambrosone - read in imputed variants in target regions

    ## read in segment file to know how many segs per chrom (made variant list per seg)
    segs <-  read.csv(file=paste(imp_dir,"/dbGaP_doctn/imputation_segments.csv",sep=""))

    vars.keep <- NA
    
    for (i in start_chr:end_chr){

      nsegs <- max(segs$segment[segs$chrom==i])

      ## read in list of imputed variant IDs in target region(s)
      for (s in 1:nsegs){
        fn.in <- paste(out_dir,"/variant_lists/vars_chr",i,".seg",s,".txt",sep="")
        tmp <- unlist(read.table(fn.in))

        vars.keep <- c(vars.keep,tmp)

      } # close segment loop
    } # close chrom loop

    ## remove type 0 not in target list from the metrics
    rm <- mets.all$type==0 & !is.element(mets.all$rs_id, vars.keep)

    ## remove from metrics files the imputed variants not in the target region(s)
    mets.all <- mets.all[!rm,]
    
    ####### Graph and summarize final imputed metrics
    if (is.element("imputed",sets)){
            ## print SNP summary
            if (summary){
              cat("Printing SNP summary by chrom\n")
              ## total study snps

              ## read in initial SNP keep list, in case type 3 SNPs were excluded from IMPUTE2 output
              kept <- read.table(file=paste(imp_dir, "/keeplists/snp.qualfilter.txt",sep=""),
                                 col.names=c("snp", "chr"), comment.char="")
              ## count by chrom
              bychr <- as.data.frame(table(kept$chr))

              ## change chrom from factor to character
              bychr$Var1 <- as.character(bychr$Var1)

              ## convert chrom to numeric code
              bychr$Var1[bychr$Var1 == "X"] <- 23
              bychr$chrom.int <- as.numeric(bychr$Var1)
              bychr <- bychr[is.element(bychr$chrom.int, start_chr:end_chr),]
              dat0 <- bychr[order(bychr$chrom.int),]

              ## total input SNPs (imp basis + study only)
              dat1 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type,c(2,3))]))
              ## total imputation basis SNPs
              dat2 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type,2)]))
              ## total imputed SNPs
              dat3 <- as.data.frame(table(mets.all$chr[is.element(mets.all$type,c(0,2))]))

              dat1$study.snps.fmkeeplist <- dat0[,2]
              dat1$imp.basis.snps <- dat2[,2]
              dat1$imp.tot.snps <- dat3[,2]
              names(dat1)[1:2] <- c("chrom", "study.snps") 
              filo <- paste(out_dir, project,"_chr",start_chr,"-",end_chr,".SNPsummary.csv",sep="")
              write.csv(dat1, file=filo, eol="\n", quote=FALSE, row.names=FALSE)
            }


            # create base name for all imputed metrics plots
            filo <- paste(out_dir, "AllImpMetrics_chr",start_chr,"-",end_chr,sep="")
            
            ## Make MAF column
            mets.all$maf <- ifelse(mets.all$exp_freq_a1 > 0.5, 1-mets.all$exp_freq_a1, mets.all$exp_freq_a1)

            ######## For Ambrosone - also make histogram of observed allele frequencies for type 2 SNPs
            cat("Plotting MAF histogram for",sum(mets.all$type==2),"type 2 SNPs\n")

            # png(paste(filo, "hist_type2MAF.png",sep="."))
            # pdf(paste(filo, "hist_type2MAF.pdf",sep="."))
            fn.out <- paste(filo, "hist_type2MAF.pdf",sep=".")
            ggplot(mets.all[mets.all$type==2,], aes(x=maf)) +
              geom_histogram(binwidth=.01, colour="black", fill="white")
            ggsave(fn.out)


            ## make data frame of all type 0 (imputed) variants -- these should be only ones in target region(s)
            mets.all.imp <- mets.all[mets.all$type==0,]

            ## ## make data frame of type 0 SNPs
            ## mets.snp.imp <- mets.all.imp[!mets.all.imp$sv.indel,]
            ## nsnps <- length(mets.snp.imp$rs_id)

            ## ## make data frame of type 0 indels and SVs
            ## mets.svi.imp <- mets.all.imp[mets.all.imp$sv.indel,]
            ## nsvi <- length(mets.svi.imp$rs_id)

            nvar <- length(mets.all.imp$rs_id)

            if(nvar < 100)
              stop("Come back when you have more imputed variants to graph, not just ",nvar,"\n")

            cat("Plots will represent a total of", nvar, "imputed variants  \n")

            ## isolate metrics column to report filtered pcts
            met <- mets.all.imp$info

            cat("\nSummaries of imputed variants passing various 'info' score thresholds:\n")

            ## loop through different thresholds
            for (t in c(0.3,0.5,0.8)){
             # cat("Fraction of SNPs with info >",t,"is",round(sum(met[!mets.all.imp$sv.indel]>t)/nsnps,4), "\n")
             # cat("Fraction of SVs and indels with info >",t,"is",round(sum(met[mets.all.imp$sv.indel]>t)/nsvi,4), "\n")
              cat("Fraction of all imputed variants with info >",t,"is",round(sum(met>t)/nvar,4), "\n")
            } # close loop on thresholds
            

            # I. Create histogram of imputation metrics
                                        # doesn't separate SNPs from SVs and indels (as these histograms aren't usually very illustrative anyway) ## DELETED SECTION AS WE'RE NOT MAKING THESE PLOTS


         # II. Create boxplots of info metric and MAF, by chrom
         ## Summarize median info metrics by chrom with box and whisker plot
         if (is.element("boxplot",plots)){

               ## loop through variant categories
               for (set in "all") {

                 # if(set=="snp") mets.curr <- mets.snp.imp
                 # if(set=="svs_indels") mets.curr <- mets.svi.imp
                 if(set=="all") mets.curr <- mets.all.imp

                    cat("making boxplots for",set,"variants\n")
                    png.fn <- paste(filo, set, "boxplot_%d.png",sep=".")
                    png(png.fn)

                    boxplot(info ~ chr, data = mets.curr, outline=FALSE, 
                    xlab = "Chromosome", ylab = "info",
                    main = NULL,las=3,col="seashell2",cex.axis=0.9)

                    boxplot(certainty ~ chr, data = mets.curr, outline=FALSE, 
                    xlab = "Chromosome", ylab = "certainty",
                    main = NULL,las=3,col="seashell2",cex.axis=0.9)

                    boxplot(maf ~ chr, data = mets.curr, outline=FALSE, 
                    xlab = "Chromosome", ylab = "imputed MAF",
                    main = NULL,las=3,col="seashell2",cex.axis=0.9)

                    dev.off()

                  } # close loop through variant types
      } # close "if" loop on whether to plot boxplots

             # III. Plot info score by MAF bin, using MAF in the study samples
             if (is.element("info_by_maf",plots)){
             
             mafs.save <- NULL

               ## loop through variant categories
               for (set in "all") {

                # if(set=="snp") mets.curr <- mets.snp.imp
                # if(set=="svs_indels") mets.curr <- mets.svi.imp
                 if(set=="all") mets.curr <- mets.all.imp

                 cat("plotting info by MAF bin for",set,"variants\n")

                 # make int.id for imputed variants
                 mets.curr$id <- 1:nrow(mets.curr)

                 attach(mets.curr)

                ####### MAF data frame loop
                mafs.list <- seq(0,0.5,length.out=51)
                mafs <- as.data.frame(mafs.list[-51]); names(mafs) <- "maf.min"
                mafs$nsnps <- mafs$mean.info <- rep(NA, times=length(mafs$maf.min))

                for (i in 1:length(mafs$maf.min)){
                  ## Count of variants per MAF bin
                  mafs$nsnps[i] <- length(id[maf > mafs.list[i] & maf <= mafs.list[i+1]])
                  ## Mean info
                  mafs$mean.info[i] <- mean(info[info!=-1 & maf > mafs.list[i] & maf <= mafs.list[i+1]])
                }
                ####### close MAF data frame loop

                detach(mets.curr)

                ## plot values at the mid point of maf bin
                xval <- mafs$maf.min+.005

                if (fmt=="png"){
                png.fn <- paste(filo, set, "info_by_studyMAF.png",sep=".")
                png(png.fn)
                }
                if (fmt=="pdf"){
                pdf.fn <- paste(filo, set,"info_by_studyMAF.pdf",sep=".")
                pdf(pdf.fn)
                }

                ## Set margins. Default is $mar 5.1 4.1 4.1 2.1 -- increase the right margin for secondary y-axis
                par(mar=c(5.1,4.1,4.1,4.1))
                plot(xval, mafs$mean.info, ylim=c(mafs$mean.info[1], 1), bty='c',
                         pch=17,xlab='MAF of imputed variants',ylab='info score (mean per MAF bin)')
                abline(h=seq(0,1,by=0.1),lty=2,col="darkgray")
                 
                par(new=TRUE)
                plot(xval, mafs$nsnps,
                         axes=FALSE, ylab='', xlab='', bty='c',pch=4, col="darkgreen")
                axis(4, tck=0.01, las=2)
                mtext("variant count", 4, col="darkgreen")
                dev.off()
                
                # save data frame
                mafs$type <- set
                mafs.save <- rbind(mafs.save,mafs)


               } # close loop on variant type

                # save MAF data frames to disk
                save(mafs.save, file=paste(filo,"info_by_MAF.RData",sep=""))

            } # close "if" loop on whether to plot info by MAF bin

    } # close "if" loop on whether to process final imp metrics

    ############ looking at imputed variants only - deleted section to Graph and summarize masked SNP metrics

} # close function definition
