#' ·············································································
#' Building SLE-Diseaseome
#' R version 4.5.1 (2025-06-13 ucrt)
#' Dec 2025
#' danieltorodominguez@gmail.com
#' ·············································································
#' Diseaseome construction

## ······································································ Step 0
## Set environment ---- 

set.seed(12345)
# memory.limit(size = 1600000000)
setwd("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/")

source("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/Code/utils.R")

check.packages(c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr",
                 "tidyr","doParallel","caret","pbapply","BiocParallel","tibble",
                 "pathMED","NbClust","ConsensusClusterPlus","SNFtool","BloodGen3Module","purrr",
                 "igraph","pheatmap","SNFtool","UpSetR","ComplexHeatmap"))

## Load Data
load(paste0(getwd(),"/RData/Datasets.rds"))
load(paste0(getwd(),"/RData/Raw_Pathway_Database.rdata"))

## Data required:
#' @DATA nested list of dataset split into Disease and Healthy (or reference)
#' matrices
#' @dbs.all list with all pathways (containing vector of genes related with)
#' @ann.info (for pathway database annotation purposes) a dataframe with
#' "annotation_id" (name of each dataset from dbs.all), "term" (description of
#' each annotation_id) and "source" (original source of each annotation)


## ······································································ Step 1
## Step 1: Increase pathway granularuty using dissectDB from pathMED ---- 

## DissectDB: split pathways into co-expressed sub-paths
#' Use only datasets with high shared genes across them (to avoid to limit the 
#' universe of genes used)

custom.db<-dissectDB(refObject = DATA[-c(4,8,12,14)], 
                     geneSets = dbs.all,
                     minPathSize = 8,
                     minSplitSize = 3,
                     maxSplits = NULL,
                     explainedVariance = 70,
                     percSharedGenes = 90) 

print(length(custom.db)) ## 47337


## Create a full table to store co-annotation information
ann.db<-as.data.frame(do.call("rbind",lapply(1:length(custom.db),function(p){
  
  nme<-gsub("\\.split.*", "", names(custom.db)[p])
  tmp<-ann.info[ann.info$annotation_id %in% nme,]
  tmp<-c(names(custom.db)[p],as.character(tmp))
  return(tmp)
})))
colnames(ann.db)<-c("TraitID","AnnotationID","term","source")
ann.db$coannotation<-NA


save.image(paste0(getwd(),"/RData/DB_step1.RData"))


## ······································································ Step 2
## Reducing pathway overlapping using Set theory (Set Packing) ---- 

# Run: 'Set Environment' section
# load("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step2.RData")

## 2.1. Remove pathways with <= 2 genes
custom.db <- custom.db[!names(custom.db) %in% 
                         names(custom.db[lengths(custom.db) <=2])]
gc()
print(length(custom.db))
## 38797

#' Get M-scores, used in setPackingFilter function, to choose between keep small
#' gene sets (gns) of the large gene set (G) each time, based on % of gained
#' information
#' 

#DATA$GSE88887$Disease<-as.matrix(DATA$GSE88887$Disease)
#DATA$GSE88887$Healthy<-as.matrix(DATA$GSE88887$Healthy)

SCORES<-mScores_createReference(refObject=DATA,
                                geneSets=custom.db,
                                cores = 16)

save.image(paste0(getwd(),"/RData/DB_step1.RData"))


## Run Set packing
time0<-Sys.time()
res.sPack<-setPackingFilter(pathDB=custom.db,
                            mRef=SCORES,
                            ann.db=ann.db,
                            coverage.thr=0.8,
                            max.combs =8, 
                            gainIndex=10,
                            minCorr=0.75)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) # 243.8337 mins

custom.db<-res.sPack$db
ann.db<-res.sPack$ann

rm(res.sPack,time0,time1)

print(length(custom.db))
## 37537

save.image(paste0(getwd(),"/RData/DB_step2.RData"))
# save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step3.RData")

## ······································································ Step 3
## Reducing pathway overlapping (Jaccard index) ---- 


## 1. Similarity 0.8 
time0<-Sys.time()
res.JI<-getNodes(listDB = custom.db,
                 ann.db=ann.db,
                 simmilarity.threshold = 0.8,
                 max.length = NULL)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) ##

custom.db<-res.JI$db
ann.db<-res.JI$ann

rm(res.JI,time0,time1)


## 2. Similarity for Small Genesets (4-5) 
time0<-Sys.time()
res.JI<-getNodes(listDB = custom.db, ann.db=ann.db,
                 simmilarity.threshold = 0.7, max.length = 4)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) ##

custom.db<-res.JI$db
ann.db<-res.JI$ann

print(length(custom.db)) ## 36853 

## Remove duplicated values in coannotation
for(i in 1:nrow(ann.db)){
  if(!is.na(ann.db[i,"coannotation"])){
    x<-unique(unlist(strsplit(ann.db[i,"coannotation"],split = ", ")))
    ann.db[i,"coannotation"]<-paste(x, collapse = ", ")
  }
}


## Filter SCORES based on current custom.db
SCORES$geneSets<-custom.db

for(i in 1:length(SCORES$mscores)){
  SCORES$mscores[[i]]<-SCORES$mscores[[i]][names(SCORES$geneSets),]
}


rm(list=setdiff(ls(),c("ann.db","custom.db","DATA","SCORES","demo.info")))


## Save R Object
save.image(paste0(getwd(),"/RData/DB_step3.RData"))


## ······································································ Step 4
## Get Disease-Relevant Gene sets (DRGs) ----

# load(paste0(getwd(),"/RData/DB_step3.RData"))
library("pathMED")

#' Step 1: "min_datasets" and "perc_samples" parameters must be optimized
#' Run 03.1_OptimizeParameters.R
#' Criteria for optimization: 
#' FDR<5%
#' LossInfo < 50 %
#' datasets > 1/3 total
#' %patietns > 10%
#' Lowest combo with same Shannon Index

#' SELECTION: min_datasets=7, perc_samples=12
## 7 Datasets, 12 % of patients

DRGs.DB<-mScores_filterPaths(MRef=SCORES,
                             min_datasets=7, perc_samples=12,
                             Pcutoff=0.05,plotMetrics = FALSE)


# Selected Genesets: 5212 genesets

full.DB<-custom.db

## Save R Object with only essential elements
rm(list=setdiff(ls(),c("full.DB","ann.db","DRGs.DB","DATA","SCORES")))
gc()

save.image(paste0(getwd(),"/RData/DB_step4.RData"))

















