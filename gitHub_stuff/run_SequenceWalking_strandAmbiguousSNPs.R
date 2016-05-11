# Function to take flanking DNA sequence of strand-ambiguous SNPs and determine TOP/BOT status
# SN 4/11/2014

# input: (1) dataframe and (2) index for column
# output: 2 columns appended to df: "myIlmnStrand" with TOP or BOT, and "myIlmnStrand.nPair" with an integer indicating the bp distance away from the SNP at which the first unambiguous pairing occured

######################

seq.walking <- function(dat, col.use, test=FALSE)
{

  require(stringr)

  # right out strings of strand ambiguous and strand unambiguous pairings
  amb.pairs <- c("T A","A T","C G","G C",
               "A A","C C","G G","T T")
  unamb.pairs <- c("A G","A C",
                "C A","C T",
                "G A","G T",
                "T G","T C")

  # if test=TRUE, run in first 100 records only
  if(test){
    cat("Testing on first 100 records only\n")
    dat.all <- dat
    dat <- dat.all[1:100,]
  }

    dat$myIlmnStrand <- NA
    dat$myIlmnStrand.nPair <- NA

    cat("Using sequence from field",names(dat)[col.use],"\n")

    ## loop over SNPs
    for (i in 1:nrow(dat)){
    # for (i in 1:100){ 

        # make uppercase version of sequence
        tmp <- toupper(dat[i,col.use])
        len <- nchar(tmp)

        # apparently the 'slash' is not always in the middle
        slash.list <- str_locate_all(pattern ='/',tmp)
        slash.pos <- slash.list[[1]][1]
        slash.chk <- substr(tmp, slash.pos, slash.pos)
        stopifnot(slash.chk=="/")

        # loop through pairs of flanking nucleotides
        # make variable that will stay 0 until we get a strand unambiguous pairing
        tk <- 0
        j <- 1
        while (tk==0) {

          # set value to add or substract
          val <- j+2
          upstr <- substr(tmp, slash.pos-val, slash.pos-val)
          downstr <- substr(tmp, slash.pos+val, slash.pos+val)

          pair <- paste(upstr,downstr)

          # decide whether to proceed - is it an ambiguous pairing?
          if(is.element(pair,amb.pairs) & !is.element(pair, unamb.pairs))
            {
            # advance to next bp
            j <- j+1
          }

          # is it an unambiguous pairing?
          if(!is.element(pair,amb.pairs) & is.element(pair, unamb.pairs))
            {
            # stop at this bp
              tk <- 1

            # assign to TOP or BOT
              if(is.element(upstr, c("A", "T"))) seq.top.bot <- "TOP"
              if(is.element(downstr, c("A", "T"))) seq.top.bot <- "BOT"
            }

          # is it something weird? (like IUPAC ambiguity code) - skip over
          if(!is.element(pair,amb.pairs) & !is.element(pair, unamb.pairs))
            warning("The pairing is wierd at sequence ",i," position",j)
        } # close loop through flanking sequence

        dat$myIlmnStrand[i] <- seq.top.bot
        dat$myIlmnStrand.nPair[i] <- j

    } # close loop over SNPS

  cat("Done!\n")

  return(dat)

} # close function definition
