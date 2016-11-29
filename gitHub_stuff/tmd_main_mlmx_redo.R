################
# MAIN ANALYSIS MLMX
################
# tmd_and MLMX
rm(list=objects())
library(OLGApipeline)
library(plyr)
sessionInfo() #OLGApipeline_0.14.1 GWASTools_1.17.8  QCpipeline_0.10.4
base.path<-"/projects/geneva/gcc-fs2/OLGA"
db <- getDb("olga_analysis")

aid.and<-640799  #tmd_and_autosomal
aid.or <-282680  #tmd_or_autosomal

db <- getDb("olga_analysis")

# since the outcome, covars, covar-matrix, sampling weights, remain same as the orginal analysis use "dbMakePrepFromAnalysis" to get these. Will have to update scan_include and config file
 
dir<-"/projects/geneva/gcc-fs2/OLGA/working_groups/dental/analysts/jaind/GWAS/analysis_prep/tmd/freeze3xup/mlmx"
dbMakePrepFromAnalysis(db, aid.and, outPath=file.path(dir,"tmd_and_redo"), stratified=FALSE)
dbMakePrepFromAnalysis(db, aid.or, outPath=file.path(dir,"tmd_or_redo"), stratified=FALSE)

# since the outcome, covars, covar-matrix, sampling weights, remain same as the orginal analysis use "dbMakePrepFromAnalysis" to get these. Will have to update scan_include and config file

# outcome
ot<-getobj(file.path(dir,"tmd_and_redo/outcome.RData")) # no changes

# covariance matrix
covm<-getobj(file.path(dir,"tmd_and_redo/covar_matrices.RData")) # no changes
#relatdness matrix addition is controlled at config file level


# covariates
cova<-getobj(file.path(dir,"tmd_and_redo/covars.RData")) # remove sex since this is sex specific analysis


# sampling weights
sw<-getobj(file.path(dir,"tmd_and_redo/weights.RData"))  # no changes

#----------------
# update scan_include
#----------------

# write scan include file
id<-getobj(file.path(dir,"tmd_and/scan_include.RData")) # scan include from the analyis were sex was inadvertently excluded
length(id) #10143
save(id, file=file.path(dir,"tmd_and_redo/scan_include.RData"))

#-----------------------
# update config
#-----------------------
cf<-read.table(file=file.path(dir,"tmd_and_redo/prep.config"),as.is=TRUE)
cf$V2[cf$V1=="analysis_description"]<-"TMDand_main effects MLM-X model"
cf$V2[cf$V1=="analysis_prefix"]<-"tmd_and_main_effects_mlmx"
cf$V2[cf$V1=="analysis_sample_description"]<-"no samples with missing values for tmd_and or eigenvectors"
a<-c("xev_set_id","17")
b<-c("xev_nums","1:2")
cf<-rbind(cf,a,b)
cf
write.table(cf, file=file.path(dir,"tmd_and_redo/prep.config"), col.names=FALSE,row.names=FALSE)

#---------
# run analysis
#---------
# prep
cd /projects/geneva/gcc-fs2/OLGA/working_groups/dental/analysts/jaind/GWAS/analysis_prep/tmd/freeze3xup/mlmx/tmd_and_redo
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/prep.py --dataset freeze3xup prep.config

#------------
# tmd_or main effects MLM-X model
#------------
# outcome
ot<-getobj(file.path(dir,"tmd_or_redo/outcome.RData")) # no changes

# covariance matrix
covm<-getobj(file.path(dir,"tmd_or_redo/covar_matrices.RData")) # no changes
#relatdness matrix addition is controlled at config file level

# covariates
cova<-getobj(file.path(dir,"tmd_or_redo/covars.RData")) # remove sex since this is sex specific analysis


# sampling weights
sw<-getobj(file.path(dir,"tmd_or_redo/weights.RData"))  # no changes


# write scan include file
id<-getobj(file.path(dir,"tmd_or/scan_include.RData")) # scan include from the analyis were sex was inadvertently excluded
length(id) #11963
save(id, file=file.path(dir,"tmd_or_redo/scan_include.RData"))

#-----------------------
# update config
#-----------------------
cf<-read.table(file=file.path(dir,"tmd_or_redo/prep.config"),as.is=TRUE)
cf$V2[cf$V1=="analysis_description"]<-"TMDor_main effects MLM-X model"
cf$V2[cf$V1=="analysis_prefix"]<-"tmd_or_main_effects_mlmx"
cf$V2[cf$V1=="analysis_sample_description"]<-"no samples with missing values for tmd_or or eigenvectors"
a<-c("xev_set_id","17")
b<-c("xev_nums","1:2")
cf<-rbind(cf,a,b)
cf
write.table(cf, file=file.path(dir,"tmd_or_redo/prep.config"), col.names=FALSE,row.names=FALSE)

#---------
# run analysis
#---------
# prep
cd /projects/geneva/gcc-fs2/OLGA/working_groups/dental/analysts/jaind/GWAS/analysis_prep/tmd/freeze3xup/mlmx/tmd_or_redo
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/prep.py --dataset freeze3xup prep.config



#--------------------------
# Check scanAnnotation for tmd_and MLM
#--------------------------
#tmd_and  MLMX scanAnnot
scanAnnot<-getobj(file.path( dir,"tmd_and_redo/scanAnnot_tmp.RData"))
dim(scanAnnot) #12803        22
head(pData(scanAnnot))
table(scanAnnot$assoc.request,scanAnnot$assoc.include)
#        FALSE  TRUE
#  FALSE  2660     0
#  TRUE      0 10143

table(scanAnnot$tmd_and[scanAnnot$assoc.include],exclude=NULL)
#   0    1 <NA>
#9374  769    0


#tmd_or  MLMX 
scanAnnot<-getobj(file.path(dir,"tmd_or_redo/scanAnnot_tmp.RData"))
dim(scanAnnot) #12803        22
head(pData(scanAnnot))
table(scanAnnot$assoc.request,scanAnnot$assoc.include)
#        FALSE  TRUE
#  FALSE   840     0
#  TRUE      0 11963

table(scanAnnot$tmd_or[scanAnnot$assoc.include],exclude=NULL) 
#   0    1 <NA>
#9374 2589    0

#-----------------
# add and submit analysis
#-----------------
## tmd_and
# add tmd_and analysis
cd /projects/geneva/gcc-fs2/OLGA/working_groups/dental/analysts/jaind/GWAS/analysis_prep/tmd/freeze3xup/mlmx/tmd_and_redo
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/prep.py --dataset freeze3xup --add prep.config
# submit assoc jobs
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/assoc.py 255422 #tmd_and main effects mlmx


## tmd_or
# add analysis
cd /projects/geneva/gcc-fs2/OLGA/working_groups/dental/analysts/jaind/GWAS/analysis_prep/tmd/freeze3xup/mlmx/tmd_or_redo
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/prep.py --dataset freeze3xup --add prep.config
# submit assoc jobs
/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/assoc.py  650686 #tmd_or  main effectsMLM_X

#-----------------
# locus zoom plot for the X chrom hit which is GWS in autosomal model
#-----------------
system(paste("/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/locus_zoom.py  --email jaind@uw.edu", 255422 ,52133577, sep=" ")) # tmd_and MLMX model, X chromsome hit
system(paste("/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/locus_zoom.py  --email jaind@uw.edu", 255422 ,46662215, sep=" ")) # tmd_and MLMX model, chr 17 hit
system(paste("/projects/geneva/gcc-fs2/OLGA/pipeline/PipelineCode/locus_zoom.py  --email jaind@uw.edu", 640799 ,46662215, sep=" ")) # tmd_and autosomal model, chr 17 hit


