#' ·············································································
#' Building SLE-Diseaseome
#' R version 4.5.1 (2025-06-13 ucrt)
#' Dec 2025
#' danieltorodominguez@gmail.com
#' ·············································································
#' Compilation of biological pathway databases / signatures


## ······································································ Step 0
## Set environment ---- 

set.seed(12345)
setwd("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/")

library("readxl")
library("stringr")
library("stringi")
library("BloodGen3Module") # Version 4.5
library("pathMED")
library("matrixStats")
library("tidyr")
library("purrr")


## ······································································ Step 1
## Collect all pathway databases ---- 

## 1. Collect pathways from databases stored in pathMED
# TMOD - Version
# Gene Ontology (BP, CC and MF) - Version
# Reactome - Version
# KEGG - Version
# Wikipathways - Version

data(genesetsData)
dbs<-genesetsData[c("tmod","go_bp","go_cc","go_mf","reactome","kegg","wikipathways")]
rm(genesetsData)

# Annotation reference object
ann.info<-do.call("rbind",lapply(1:length(dbs),function(i){
    paths<-names(dbs[[i]])
    tmp<-pathMED:::ann_info
    tmp<-tmp[tmp$annotation_id %in% paths,]
    tmp$source<-rep(names(dbs)[i],nrow(tmp))
    return(tmp)
}))


## 2. Cell type signatures (from Single Cell)

#' Get single cell Genesets FROM xCellDB
## Download xCell database from https://pmc.ncbi.nlm.nih.gov/articles/PMC5688663/
# Aditional File 3 (489 Inicial Signatures)
data <- read_excel("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/Datasets/xCellDB.xlsx")
data<-as.data.frame(data)
data$Celltype_Source_ID<-gsub("\\+","pos",data$Celltype_Source_ID)

uniqPaths<-data$Celltype_Source_ID
uniqPaths<-unique(stringi::stri_replace_all_regex(uniqPaths,
                                                  pattern = c("_1", "_2", "_3","_HPCA","_ENCODE","_FANTOM","_IRIS","_BLUEPRINT","_NOVERSHTERN"), 
                                                  replacement = c("","","","","","","","",""),
                                                  vectorize = FALSE))

xcelldb<-lapply(1:length(uniqPaths),function(i){
  tmp<-data[str_detect(data$Celltype_Source_ID,pattern = uniqPaths[i]),]
  tmp<-tmp[,-c(1:2)]
  tmp<-unlist(tmp)
  tmp<-tmp[!is.na(tmp)]
  return(unique(tmp))
})
names(xcelldb)<-gsub(" ","_",uniqPaths)
names(xcelldb)<-gsub("-","_",names(xcelldb))
## 64 Signatures after remove redundancies


## Add 10 additional Cell Signatures from bibliography curated by experts
# https://doi.org/10.3389/fimmu.2021.602539 - PMC8012727
# DOI:10.1038/s41590-024-01883-0 - PMC11291291

bibl<-list("B_cells_DN4"=c("HOPX","PDE4D","IGHE","SELL"),
           "B_cells_DN3"= c("RHOB","VIM","CTSH"),
           "B_cells_DN2"= c("EMP3","CIB1","PSAP","CD72","DAPP1","HCK","ZEB2","RHOB","TNFRSF1B","FCRL3","FCRL5","FGR","MPP6"),
           "B_cells_DN1" = c("TAGLN2","IGHA2","JCHAIN","IGHA1","S100A10"),
           "B_cells_trans"=c("VPREB3","IGHD","IIGLL5","TCL1A"),
           "APCs" = c("IFI30","TNFAIP2","CLEC7A"),
           "NK_cells_NKT" = c("CD3E","CD3D","CD3G","HLA-DPB1","HLA-DRB1","HLA-DQB1"),
           "NK_cells_CD56bright" = c("GPR183","ILR7","LTB","GZMK","CD62L","CCR7","CD2","KLRC1"),
           "NK_cells_CD56dim_CD16pos" = c("FCGR3A","KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2","KIR3DL3","KIR2DL4"),
           "NK_cells_CD56neg_CD16pos_CD7pos" = c("CCL3","CCL4","CCL5"))


nonBlood<-c("Astrocytes","Hepatocytes","Adipocytes","aDC","Chondrocytes",
            "Endothelial","Epithelial","Keratinocytes","Melanocytes",
            "Mesangial","MSC","Neurons","Osteoblast","Preadipocytes","muscle")

# Filter non-blood cells
xcelldb<-xcelldb[!sapply(names(xcelldb), function(element, pattern) {
  sapply(pattern, function(patr) grepl(patr, element)) %>% any()
}, nonBlood)]

ann.info<-rbind(ann.info,data.frame("annotation_id"=names(xcelldb),
                                    "term"=names(xcelldb),
                                    "source"="xcell"))

xcelldb<-append(xcelldb,bibl)
## 56 Final Blood Signatures from Single Cell (xCell + 2 articles)

ann.info<-rbind(ann.info,data.frame("annotation_id"=names(bibl),
                                    "term"=names(bibl),
                                    "source"=c(rep("PMC8012727",6),rep("PMC11291291",4))))

## Save Object
saveRDS(xcelldb,"C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/xcellDB.rds")

dbs[["xcell"]]<-xcelldb
rm(data,uniqPaths,nonBlood,xcelldb,bibl)


## 3. Add SLE-relevant signatures from the Bibliography
# DOI: 10.1016/j.xcrm.2024.101569 - PMC11148857

ifnSig<-read.csv("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/Datasets/IFN_Signature.csv",
                 sep=";")[,1:2]
ifnSig <- as.data.frame(ifnSig %>% separate_rows(Gene, sep = " /// "))

knowlge<-list("IFN_I" = as.character(ifnSig[ifnSig$Signature=="IFI_I","Gene"]),
              "IFN_II" = as.character(ifnSig[ifnSig$Signature=="IFI_II","Gene"]),
              "IFN_II_III" = as.character(ifnSig[ifnSig$Signature=="IFI_II_III","Gene"]),
              "IFN_I_II_III" = as.character(ifnSig[ifnSig$Signature=="IFI_I_II_III","Gene"]),
              "T_cell_exhaustion" = c("CTLA4", "IL7R", "LAG3", "PDCD1", "ABCE1"))

ann.info<-rbind(ann.info,data.frame("annotation_id"=names(knowlge),
                                    "term"=names(knowlge),
                                    "source"="PMC11148857"))

dbs[["knowlge"]]<-knowlge
rm(ifnSig,knowlge)


## 4. Add signatures from BloodGen3Module - Version 4.5
annb3m<-as.data.frame(BloodGen3Module:::Module_listGen3)

b3m<-lapply(unique(annb3m$Module),function(term){
  tmp<-annb3m[annb3m$Module==term,]
  return(as.character(tmp$Gene))
})
names(b3m)<-unique(annb3m$Module)
dbs[["B3M"]]<-b3m

ann.b3m<-do.call("rbind",lapply(1:length(b3m),function(i){
    paths<-names(b3m)[i]
    tmp<-annb3m[annb3m$Module %in% paths,c("Module","Function")]
    tmp<-tmp[!duplicated(tmp),]
    colnames(tmp)<-c("annotation_id","term")
    return(tmp)
}))
ann.b3m$source<-"B3M"

ann.info<-rbind(ann.info,ann.b3m)
print(table(ann.info$source)) # 23676

## merge all databases
dbs.all<-flatten(dbs)

rm(b3m,ann.b3m,annb3m,dbs)

save.image("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/RData/Raw_Pathway_Database.rdata")




