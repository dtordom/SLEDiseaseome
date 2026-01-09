#' ·············································································
#' Building SLE-Diseaseome
#' R version 4.5.1 (2025-06-13 ucrt)
#' Dec 2025
#' danieltorodominguez@gmail.com
#' ·············································································
#' Compilation of Gene-Expression Datasets from SLE and NHV samples


## ······································································ Step 0
## Set environment ---- 

set.seed(12345)
# options(timeout = 1e5)
setwd("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/")

source("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/Code/utils.R")

check.packages(c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr",
                 "tidyr","doParallel","caret","pbapply","BiocParallel","tibble",
                 "GEOquery","data.table","ggplot2","GEOmetadb","readxl"))

## Save ensembl to annotate genes
ensembl = useEnsembl(biomart='ensembl', dataset="hsapiens_gene_ensembl",
                     mirror = "useast") # useast, www,asia
# Save Object to avoid posible future server errors
# saveRDS(ensembl,"C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/Ensembl.rds")

## List to store the datasets
DATA<-list()
demo.info<-list()
#' demo.info: 
#' Disease Activity, Gender, Age, Race, Geographic region, Tissue


## ······································································ Step 1
## Load data from PRECISESADs ---- 
#' 2 sets of patients and 1 of healthy controls, without batch effect)
#' Available by request to the authors: https://pubmed.ncbi.nlm.nih.gov/33497037/

Count.precisesads<-read.table(paste0(getwd(),"/Datasets/Count.PRECISESADS.csv"),
                      header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads<-type.convert(x = Count.precisesads,as.is=F)

## Filter non-expressed genes
Count.precisesads<-Count.precisesads[rownames(Count.precisesads)[rowCounts(Count.precisesads>=10)>=10],]

## Annotate to Gene symbol
Count.precisesads<-annotateGenes(data=Count.precisesads,
                                 toGenes='external_gene_name',
                                 fromGenes='ensembl_gene_id',
                                 method = "median",ensembl = ensembl)

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads)
data<-log2(data+1)
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads1<-data
rm(Count.precisesads,data)


Count.precisesads2<-read.table(paste0(getwd(),"/Datasets/Count.PRECISESADS2.csv"),
                               header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads2<-type.convert(x = Count.precisesads2,as.is=F)

## Filter non-expressed genes
Count.precisesads2<-Count.precisesads2[rownames(Count.precisesads2)[rowCounts(Count.precisesads2>=10)>=10],]

## Annotate to Gene symbol
Count.precisesads2<-annotateGenes(data=Count.precisesads2,
                                  toGenes='external_gene_name',
                                  fromGenes='ensembl_gene_id',
                                  method = "median",ensembl = ensembl)

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads2)
data<-log2(data+1)
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads2<-data
rm(Count.precisesads2,data)


## Select SLE and Healthy samples 
clin1<-read.csv(file=paste0(getwd(),"/Datasets/Metadata.PRECISESADS.csv"),
                header=T,sep=";",row.names=1,dec=",")
clin1<-clin1[colnames(precisesads1),c("Diagnosis","Age","Gender","Center")]

SLE<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="SLE",T,F),])]
HC.prec<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="CTRL",T,F),])]  

## Remove non-variable genes in HC or SLE independently
#' To avoid technical issues in specific groups
## NHV
nonVar.genes<-summ.var(data = HC.prec,freqCut = 70/30,
                       uniqueCut = 30,minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC.prec<-HC.prec[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30,minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 


clin2<-read.csv(file= paste0(getwd(),"/Datasets/Metadata.PRECISESADS2.csv"),
                header=T,sep=";",row.names=1,dec=",")
clin2<-clin2[colnames(precisesads2),c("Diagnosis","Age","Gender","Center")]

SLE2<-precisesads2[,rownames(clin2[ifelse(clin2$Diagnosis=="SLE",T,F),])]


## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = SLE2,freqCut = 70/30,
                       uniqueCut = 30,minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE2<-SLE2[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## Select common genes
genes<-list(rownames(HC.prec),rownames(SLE),rownames(SLE2))
genes<-Reduce(intersect,genes)

SLE.prec<-cbind(SLE[genes,],SLE2[genes,])
HC.prec<-HC.prec[genes,]

DATA[["precisesads"]]<-list("Disease"=SLE.prec,"Healthy"=HC.prec)

rm(clin1,clin2,SLE,genes,SLE2,precisesads1,precisesads2,HC.prec,SLE.prec,nonVar.genes)


## Load Demographical data (PRECISESADS)
clin<-readRDS(paste0(getwd(),"/Datasets/Precisesads_SLEDAI.rds"))

clin<-clin[clin$OMICID %in% intersect(colnames(DATA$precisesads$Disease),clin$OMICID),]
rownames(clin)<-clin$OMICID

clin$GeographicRegion<-"EU"

demo.info[["precisesads"]]<-clin
rm(clin)


## ······································································ Step 2
## Load data from Pascual et. al (GSE65391) ----
## gene expression for GSE65391 can be also downloaded from GEO using getGEO function

data<-read.table(file=paste0(getwd(),"/Datasets/GSE65391_expressionData.txt"),
                 sep="\t",header=T,row.names = 1)
clin<-data.frame(t(read.csv(file=paste0(getwd(),"/Datasets/Pascual_allClin.csv"),
                            sep="\t",row.names = "Sample_ID")))
data<-data[,rownames(clin)];

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)

## Separate Healthy and SLE samples
clin$state<-ifelse(grepl("BAY", clin$Patient),"Healthy","SLE")
HC<-data[,clin$state=="Healthy"]
SLE<-data[,clin$state=="SLE"]

## Remove non-variable genes in HC or SLE independently
##NHV
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 10, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 10, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 


genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["pascual"]]<-list("Disease"=SLE,"Healthy"=HC)

age<-getGEO("GSE65391", GSEMatrix =TRUE,
              destdir=paste0(getwd(),"/RData"))
age<-phenoData(age[[1]])
age<-pData(age)
age<-data.frame("SampleID"=rownames(age),
                "Age"=as.numeric(age$'age:ch1'))
rownames(age)<-age$SampleID
age<-age[rownames(clin),]

## Get mean SLEDAI per patient
da<-unlist(lapply(clin$Patient,function(p){
  return(mean(as.numeric(clin[clin$Patient == p,"SLEDAI"]),na.rm=T))
}))

clin<-data.frame("SampleID"=rownames(clin),
                 "PatientID"=clin$Patient,
                 "SLEDAI"=da,
                 "Age"=age$Age,
                 "Gender"=clin$Gender,
                 "Race"=clin$Race,
                 "GeographicRegion"="USA",
                 "Tissue"="WholeBlood")
rownames(clin)<-clin$SampleID
clin<-clin[colnames(SLE),]
clin<-clin[!duplicated(clin$PatientID),]

demo.info[["pascual"]]<-clin

rm(clin,data,HC,SLE,nonVar.genes,genes,da,age)

# save.image("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/tmpRData.rdata")



#### Jonhs Hopkins cohort ----------------------------------------------- Step 3
## All clinical data can be requested to the original authors (GSE45291)

data<-read.csv(file=paste0(getwd(),"/Datasets/PetriALL.txt"),
               sep="\t",row.names = "GeneSymbol")

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

## Separate Healthy and SLE samples
clin<-read.csv(file=paste0(getwd(),"/Datasets/Metadata.petri.csv"),
               sep=";",row.names = "GZ_Filenames")

HC<-data[,rownames(clin)[clin$Diagnosis=="Healthy"]]
SLE<-data[,rownames(clin)[clin$Diagnosis=="SLE"]]

## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["petri"]]<-list("Disease"=SLE,"Healthy"=HC)

## demographical information
clin1<-read.csv(paste0(getwd(),"/Datasets/histol.petri.csv"),sep="\t")
clin2<-read.csv(paste0(getwd(),"/Datasets/da.petri.csv"),sep="\t")[,c("CohortID","SLEDAI")]

## Get mean SLEDAI per patient
da<-unlist(lapply(clin1$Cohort.ID,function(p){
  return(mean(as.numeric(clin2[clin2$CohortID == p,"SLEDAI"]),na.rm=T))
}))

clin<-data.frame("SampleID"=clin1$Cohort.ID,
                 "PatientID"=clin1$Cohort.ID,
                 "SLEDAI"=da,
                 "Age"=as.numeric(clin1$Age.Onset),
                 "Gender"=clin1$Sex.UC,
                 "Race"=clin1$Race.UC,
                 "GeographicRegion"="USA",
                 "Tissue"="PeripheralBlood")

demo.info[["petri"]]<-clin

rm(clin,clin1,clin2,data,HC,SLE,nonVar.genes,genes,da)

# save.image("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/tmpRData.rdata")


#### Datasets from ADEX ------------------------------------------------- Step 4
## https://adex.genyo.es/

clin<-read.csv(file=paste0(getwd(),"/Datasets/metadataAdex.tsv"),sep="\t")
rownames(clin)<-clin$Sample
# Tissue, Gender, Age and Race

# "GSE45291_SLE.tsv"
datasets<-c("GSE24706.tsv","GSE50772.tsv","GSE61635.tsv",
            "GSE72509.tsv","GSE82221_GPL10558.tsv","GSE108497.tsv",
            "GSE110169_SLE.tsv","GSE110174.tsv")

## Transcriptome data
for(d in 1:length(datasets)){
  dat<-datasets[d]
  datName<-gsub(".tsv","",dat)
  
  data<-read.csv(file=paste0(getwd(),"/Datasets/",dat),
                 sep="\t",row.names = "gene")
  
  ## Preprocessing (Log normalization and filtering non variable genes)
  data<-norm.log(data)
  
  ## Separate Healthy and SLE samples
  tmp<-clin[colnames(data),c("GSE","Condition")]
  HC<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
  SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]
  
  ## Remove non-variable genes in HC or SLE independently
  nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                         uniqueCut = 30, 
                         minVar = 0.05)
  table(nonVar.genes$nzv)
  table(nonVar.genes$thrSD)
  plot(density(nonVar.genes$sd))
  print(min(nonVar.genes$sd))
  HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 
  
  ## SLE
  nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                         uniqueCut = 30, 
                         minVar = 0.05)
  table(nonVar.genes$nzv)
  table(nonVar.genes$thrSD)
  plot(density(nonVar.genes$sd))
  print(min(nonVar.genes$sd))
  SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 
  
  
  genes<-intersect(rownames(HC),rownames(SLE))
  HC<-HC[genes,]
  SLE<-SLE[genes,]
  
  DATA[[datName]]<-list("Disease"=SLE,"Healthy"=HC)
}


## Demographical information for each dataset

## GSE24706
tmpClin<-clin[clin$GSE=="GSE24706",]
demo.info[["GSE24706"]]<-data.frame("SampleID"=tmpClin$Sample,
                                "PatientID"=tmpClin$Sample,
                                "SLEDAI"=NA,"Age"=NA,
                                "Gender"=NA,"Race"=NA,
                                "GeographicRegion"="USA",
                                "Tissue"="PBMC")

## GSE50772
tmpClin<-clin[clin$GSE=="GSE50772",]
demo.info[["GSE50772"]]<-data.frame("SampleID"=tmpClin$Sample,
                                    "PatientID"=tmpClin$Sample,
                                    "SLEDAI"=NA,"Age"=NA,
                                    "Gender"=NA,"Race"=NA,
                                    "GeographicRegion"="USA",
                                    "Tissue"="PBMC")

## GSE61635
tmpClin<-clin[clin$GSE=="GSE61635",]
demo.info[["GSE61635"]]<-data.frame("SampleID"=tmpClin$Sample,
                                    "PatientID"=tmpClin$Sample,
                                    "SLEDAI"=NA,"Age"=NA,
                                    "Gender"=NA,"Race"=NA,
                                    "GeographicRegion"=NA,
                                    "Tissue"="WholeBlood")

## GSE72509
tmpClin<-clin[clin$GSE=="GSE72509",]
demo.info[["GSE72509"]]<-data.frame("SampleID"=tmpClin$Sample,
                                    "PatientID"=tmpClin$Sample,
                                    "SLEDAI"=NA,"Age"=NA,
                                    "Gender"=NA,"Race"=NA,
                                    "GeographicRegion"="USA",
                                    "Tissue"="WholeBlood")

## GSE82221_GPL10558
tmpClin<-clin[clin$GSE=="GSE82221",]

tmpClin2<-getGEO("GSE82221", GSEMatrix =TRUE,
            destdir=paste0(getwd(),"/RData"))
tmpClin2<-phenoData(tmpClin2[[1]])
tmpClin2<-pData(tmpClin2)
tmpClin2<-tmpClin2[gsub("_.*","",tmpClin$Sample),c("age:ch1","Sex:ch1")]
tmpClin2$`age:ch1`<-as.numeric(as.character(gsub("y","",tmpClin2$`age:ch1`)))

demo.info[["GSE82221_GPL10558"]]<-data.frame("SampleID"=tmpClin$Sample,
                                    "PatientID"=tmpClin$Sample,
                                    "SLEDAI"=NA,
                                    "Age"=tmpClin2$`age:ch1`,
                                    "Gender"=tmpClin2$`Sex:ch1`,
                                    "Race"="Asian",
                                    "GeographicRegion"="China",
                                    "Tissue"="PBMC")

## GSE108497
tmpClin<-clin[clin$GSE=="GSE108497",]

tmpClin2<-getGEO("GSE108497", GSEMatrix =TRUE,
                 destdir=paste0(getwd(),"/RData"))
tmpClin2<-phenoData(tmpClin2[[1]])
tmpClin2<-pData(tmpClin2)

tmpClin2<-tmpClin2[tmpClin$Sample,]

demo.info[["GSE108497"]]<-data.frame("SampleID"=tmpClin$Sample,
                                    "PatientID"=tmpClin$Sample,
                                    "SLEDAI"=NA,
                                    "Age"=as.numeric(tmpClin2$`age:ch1`),
                                    "Gender"=tmpClin$Gender,
                                    "Race"=tmpClin2$`race:ch1`,
                                    "GeographicRegion"="USA",
                                    "Tissue"="WholeBlood")

## GSE110169_SLE
tmpClin<-clin[clin$GSE=="GSE110169",]

demo.info[["GSE110169_SLE"]]<-data.frame("SampleID"=tmpClin$Sample,
                                     "PatientID"=tmpClin$Sample,
                                     "SLEDAI"=NA,"Age"=NA,
                                     "Gender"=tmpClin$Gender,
                                     "Race"=NA,
                                     "GeographicRegion"="USA",
                                     "Tissue"="Peripheralblood")


## GSE110174
tmpClin<-clin[clin$GSE=="GSE110174",]

demo.info[["GSE110174"]]<-data.frame("SampleID"=tmpClin$Sample,
                                     "PatientID"=tmpClin$Sample,
                                     "SLEDAI"=NA,"Age"=NA,
                                     "Gender"=NA,"Race"=NA,
                                     "GeographicRegion"="USA",
                                     "Tissue"="Peripheralblood")  

rm(clin,data,HC,nonVar.genes,SLE,tmp,tmpClin,tmpClin2,d,dat,datasets,datName,
   genes)

# save.image("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/tmpRData.rdata")


#### Query NCBI GEO ----------------------------------------------------- Step 5
## GSE211700, GSE22098, GSE88887

##  GSE211700 ·············
## Matrix downloaded from NCBI GEO
data <- as.data.frame(fread(paste0(getwd(),"/Datasets/GSE211700_Transcript_FPKM.txt")))

genome<-data[,c("Transcript_id","Official_Symbol")]
rownames(data)<-data$Transcript_id
data<-data[,-c(1:4)]

data<-norm.log(data+1)

genome$Official_Symbol<-gsub("\\..*","",genome$Official_Symbol)
genome[genome=="--"]<-NA

genome <- genome %>% `colnames<-`(c("fromGenes","toGenes")) %>% 
  replace(.=="",NA) %>% drop_na() %>% filter(fromGenes %in% as.character(rownames(data)))

data = data[genome$fromGenes,]
finalGenes = unique(genome$toGenes)

temp = as.data.frame(do.call("rbind", mclapply(finalGenes, mergeExpression,
                                               genome = genome,
                                               expressionMatrix = data,
                                               method = "median",
                                               mc.cores = ifelse(as.character(Sys.info()["sysname"])=="Windows",1,detectCores()))))
rownames(temp) = finalGenes
colnames(temp) = colnames(data)
data<-temp

SLE<-data[,-c(grep("CTRL",colnames(data)))]
HC<-data[,grep("CTRL",colnames(data))]

## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 


genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["GSE211700"]]<-list("Disease"=SLE,"Healthy"=HC)

## Demographical info
clin<-getGEO("GSE211700", GSEMatrix =TRUE, destdir=paste0(getwd(),"/RData"))
clin<-phenoData(clin[[1]])
clin<-pData(clin)

demo.info[["GSE211700"]]<-data.frame("SampleID"=clin$geo_accession,
                                     "PatientID"=clin$geo_accession,
                                     "SLEDAI"=NA,"Age"=NA,"Gender"=NA,
                                     "Race"="Asian",
                                     "GeographicRegion"="China",
                                     "Tissue"="PBMC")  

rm(clin,temp,finalGenes,genome,nonVar.genes,SLE,HC,data,genes)



## GSE22098 ·············
gset = getGEO("GSE22098", GSEMatrix =TRUE, destdir=paste0(getwd(),"/RData"))
clin<-phenoData(gset[[1]])
clin<-pData(clin)
clin<-clin[,c("title","geo_accession","healthy control:ch1","illness:ch1",
              "gender:ch1","ethnicity:ch1","age:ch1")]
gset<-gset[[1]]
data<-exprs(gset)

data<-norm.log(data + abs(min(data))+1)

genome<-read.csv(paste0(getwd(),"/Datasets/GPL6947-13512.txt"),
                 sep="\t",skip=30)[,c("ID","Symbol")]

genome <- genome %>% `colnames<-`(c("fromGenes","toGenes")) %>% 
  replace(.=="",NA) %>% drop_na() %>% filter(fromGenes %in% as.character(rownames(data)))

data = data[genome$fromGenes,]
finalGenes = unique(genome$toGenes)

temp = as.data.frame(do.call("rbind", mclapply(finalGenes, mergeExpression,
                                               genome = genome,
                                               expressionMatrix = data,
                                               method = "median",
                                               mc.cores = ifelse(as.character(Sys.info()["sysname"])=="Windows",1,detectCores()))))
rownames(temp) = finalGenes
colnames(temp) = colnames(data)
data<-temp

colnames(clin)<-c("title","gsm","healthy","disease","gender","race","age")

## Pediatric set
hc<- clin %>% filter(healthy == "pSLE")
HC<-data[,hc$gsm]
sle<- clin %>% filter(disease == "PSLE")
SLE<-data[,sle$gsm]

## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 


genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["GSE22098p"]]<-list("Disease"=SLE,"Healthy"=HC)

demo.info[["GSE22098p"]]<-data.frame("SampleID"=clin[colnames(SLE),"gsm"],
                                     "PatientID"=clin[colnames(SLE),"gsm"],
                                     "SLEDAI"=NA,
                                     "Age"=as.numeric(clin[colnames(SLE),"age"]),
                                     "Gender"=clin[colnames(SLE),"gender"],
                                     "Race"=clin[colnames(SLE),"race"],
                                     "GeographicRegion"="USA",
                                     "Tissue"="WholeBlood")  


## Adult set
hc<- clin %>% filter(healthy == "ASLE")
HC<-data[,hc$gsm]
sle<- clin %>% filter(disease == "ASLE")
SLE<-data[,sle$gsm]

## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["GSE22098a"]]<-list("Disease"=SLE,"Healthy"=HC)

demo.info[["GSE22098a"]]<-data.frame("SampleID"=clin[colnames(SLE),"gsm"],
                                     "PatientID"=clin[colnames(SLE),"gsm"],
                                     "SLEDAI"=NA,
                                     "Age"=as.numeric(clin[colnames(SLE),"age"]),
                                     "Gender"=clin[colnames(SLE),"gender"],
                                     "Race"=clin[colnames(SLE),"race"],
                                     "GeographicRegion"="USA",
                                     "Tissue"="WholeBlood")  

rm(data,clin,genome,gset,hc,HC,SLE,sle,temp,finalGenes,nonVar.genes,genes)


#### Query NCBI GEO ----------------------------------------------------- Step 6
## GSE88887 - Tabalumab (large dataset)

load(paste0(getwd(),"/Datasets/Tabalumab.RData"))
## ONLY Baseline samples, or samples at week16 and 52 but treated with PBO, were selected

data<-norm.log(Raw.tab)

## Annotation to gene symbol
# [HTA-2_0] Affymetrix Human Transcriptome Array 2.0 [transcript (gene) version]
rownames(data)<-gsub("(.hg).*","\\1",rownames(data))

data<-annotateGenes(data=data,
                    toGenes='external_gene_name',
                    fromGenes='affy_hta_2_0',
                    method = "median",ensembl=ensembl)

HC<-data[,rownames(clin.tab[clin.tab$State=="Normal",])]
SLE<-data[,rownames(clin.tab[clin.tab$State=="SLE",])]

## Remove non-variable genes in HC or SLE independently
nonVar.genes<-summ.var(data = HC,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
HC<-HC[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 

## SLE
nonVar.genes<-summ.var(data = SLE,freqCut = 70/30,
                       uniqueCut = 30, 
                       minVar = 0.05)
table(nonVar.genes$nzv)
table(nonVar.genes$thrSD)
plot(density(nonVar.genes$sd))
print(min(nonVar.genes$sd))
SLE<-SLE[!(nonVar.genes$nzv | nonVar.genes$thrSD),] 


genes<-intersect(rownames(HC),rownames(SLE))
HC<-HC[genes,]
SLE<-SLE[genes,]

DATA[["GSE88887"]]<-list("Disease"=SLE,"Healthy"=HC)

## Get demographical info

gset = getGEO("GSE88887", GSEMatrix =TRUE, destdir=paste0(getwd(),"/RData"))
clin<-phenoData(gset[[1]])
clin<-pData(clin)

clin<-clin[rownames(clin.tab),]
clin<-clin[clin$'time:ch1'=="baseline",c("race:ch1","Sex:ch1","sledai_at_baseline:ch1",
                                         "subject_id:ch1","region:ch1","age_at_baseline:ch1")]

clin[clin=="--"]<-NA


demo.info[["GSE88887"]]<-data.frame("SampleID"=rownames(clin),
                                     "PatientID"=clin$`subject_id:ch1`,
                                     "SLEDAI"=as.numeric(clin$`sledai_at_baseline:ch1`),
                                     "Age"=as.numeric(clin$`age_at_baseline:ch1`),
                                     "Gender"=clin$`Sex:ch1`,
                                     "Race"=clin$`race:ch1`,
                                     "GeographicRegion"="USA",
                                     "Tissue"="WholeBlood")  

 rm(list=setdiff(ls(),c("DATA","demo.info")))

save.image(paste0(getwd(),"/RData/Datasets.rds"))



#### Visualization ------------------------------------------------------ Step 7
## 

library("ggplot2")
library("ggpubr")

load(paste0(getwd(),"/RData/Datasets.rds"))

## Age
lapply(1:length(demo.info),function(i){
  tmp<-demo.info[[i]]
  tmp<-tmp$Age
  mn<-round(mean(tmp,na.rm=T),digits = 2)
  sd<-round(sd(tmp,na.rm=T),digits = 2)
  return(paste0(mn," (±",sd,")"))
})

## SLEDAI
lapply(1:length(demo.info),function(i){
  tmp<-demo.info[[i]]
  tmp<-tmp$SLEDAI
  mn<-round(mean(tmp,na.rm=T),digits = 2)
  sd<-round(sd(tmp,na.rm=T),digits = 2)
  return(paste0(mn," (±",sd,")"))
})


## Gender
lapply(1:length(demo.info),function(i){
  tmp<-demo.info[[i]]
  if(i!=3){
    tmp<-tmp$Gender[tmp[,1] %in% colnames(DATA[[i]]$Disease)]
  }else{
    tmp<-tmp$Gender
  }
  
  tmp<-ifelse(tmp=="Female","F",tmp)
  tmp<-ifelse(tmp=="Male","M",tmp)
  
  if(!all(is.na(tmp))){
    
    return(paste0(table(tmp)["F"]," / ",table(tmp)["M"]))
    
  }else{
    return(NA)
  }
})

## Region
lapply(1:length(demo.info),function(i){
  tmp<-demo.info[[i]]
  tmp<-tmp$GeographicRegion
  return(unique(tmp))
})


## Race
lapply(1:length(demo.info),function(i){
  tmp<-demo.info[[i]]
  return(table(tmp$Race))
})


theme_set(theme_bw() + 
            theme(legend.position = "bottom",
                  legend.title = element_text(size=8, face="bold"),
                  legend.text = element_text(size=8),
                  axis.text.x=element_text(size=8, color = "black"),
                  axis.text.y = element_text(size=8, color = "black"),
                  axis.title = element_text(size=10, color = "black", face = "bold"), 
                  plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
                  plot.subtitle = element_text(size = 8,hjust = 0.5),
                  strip.text = element_text(size=6), 
                  strip.background = element_blank()))


m<- as.data.frame(matrix(data=0,nrow=length(DATA),ncol=3))
colnames(m)<-c("study","healthy","sle")
rownames(m)<-names(DATA)
m$study<-names(DATA)

for (i in 1:length(DATA)) {
  m$sle[i]<-ncol(DATA[[i]]$Disease)
  m$healthy[i]<-ncol(DATA[[i]]$Healthy)
}

## Longitudinal patients: pascual, petri, GSE88887, GSE108497
m$patients<-m$sle
m$patients[2]<-158
m$patients[3]<-301
m$patients[9]<-92
m$patients[15]<-1760


m$platform<-c("RNASeq","Microarray","Microarray","Microarray",
              "Microarray","Microarray","RNASeq","Microarray","Microarray",
              "Microarray","Microarray","RNASeq","Microarray","Microarray",
              "Microarray")
m$age<-c("adult","pediatric","adult","adult","adult","adult","adult",
         "adult","adult","adult","adult","adult","pediatric","adult","adult")

m$tissue<-c("WB","PBMC","WB","WB","PBMC","WB","WB",
            "PBMC","WB","WB","WB","PBMC","WB","WB","WB")


## SLE vs NHV
m.1<- data.frame("group"=c("SLE","NHV"),
                 "value"=c(sum(m$sle),sum(m$healthy)),
                 "text"=c(paste0(sum(m$sle)," (",sum(m$patients),")"),
                          paste0(sum(m$healthy)," (",sum(m$healthy),")")))

m.1 <- m.1 %>% 
  mutate(prop = value / sum(value) * 100,
         ypos = cumsum(prop) - 0.5 * prop)

p1<-ggplot(m.1, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  labs(title = "Sample size") +
  theme_void() + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold")) +
  scale_fill_manual(values=c("SLE"="lightskyblue3",
                             "NHV"= "#CCCCCC")) +
  geom_text(aes(label = text, y = ypos), color = "black", size = 3)


## Platform
m.2<-m
m.2$total<-apply(m.2[,c("sle","healthy")],1,sum)

m.2<- as.data.frame(m.2 %>% 
                      group_by(platform) %>%
                      summarise(val = sum(total)))

m.2 <- m.2 %>% mutate(prop = val / sum(val) * 100,
                      ypos = cumsum(prop) - 0.5 * prop)
m.2$platform<-factor(m.2$platform,levels=c("RNASeq","Microarray"))

p2<-ggplot(m.2, aes(x = "", y = prop, fill = platform)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  labs(title = "Platform") +
  theme_void() + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold")) +
  scale_fill_manual(values=c("RNASeq"="#CCCCCC",
                             "Microarray"= "lightskyblue3")) +
  geom_text(aes(label = val, y = ypos), color = "black", size = 3)


## Age
m.2<-m
m.2$total<-apply(m.2[,c("sle","healthy")],1,sum)

m.2<- as.data.frame(m.2 %>% 
                      group_by(age) %>%
                      summarise(val = sum(total)))

m.2 <- m.2 %>% mutate(prop = val / sum(val) * 100,
                      ypos = cumsum(prop) - 0.5 * prop)
m.2$age<-factor(m.2$age,levels=c("pediatric","adult"))

p3<-ggplot(m.2, aes(x = "", y = prop, fill = age)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  labs(title = "Age") +
  theme_void() + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold")) +
  scale_fill_manual(values=c("pediatric"="#CCCCCC",
                             "adult"= "lightskyblue3")) +
  geom_text(aes(label = val, y = ypos), color = "black", size = 3)


## tissue
m.2<-m
m.2$total<-apply(m.2[,c("sle","healthy")],1,sum)

m.2<- as.data.frame(m.2 %>% 
                      group_by(tissue) %>%
                      summarise(val = sum(total)))

m.2 <- m.2 %>% mutate(prop = val / sum(val) * 100,
                      ypos = cumsum(prop) - 0.5 * prop)
m.2$tissue<-factor(m.2$tissue,levels=c("WB","PBMC"))

p4<-ggplot(m.2, aes(x = "", y = prop, fill = tissue)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  labs(title = "Tissue") +
  theme_void() + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold")) +
  scale_fill_manual(values=c("WB"="lightskyblue3",
                             "PBMC"= "#CCCCCC")) +
  geom_text(aes(label = val, y = ypos), color = "black", size = 3)


ggarrange(p1,p2,p3,p4,ncol=4,nrow = 1)






