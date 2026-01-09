#' ·············································································
#' Building SLE-Diseaseome
#' R version 4.5.1 (2025-06-13 ucrt)
#' Dec 2025
#' danieltorodominguez@gmail.com
#' ·············································································
#' Optimize parameters for mScores_filterPaths function (pathMED)
#' @min_datasets number of datasets that each pathway must meet the perc_samples
#' threshold
#' @perc_samples minimun percentage of samples in a dataset in which a pathway
#' must be significant

## ······································································ Step 0
## Set environment ---- 

set.seed(12345)
# memory.limit(size = 1600000000)
setwd("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/")

source("C:/Users/danie/Desktop/WORK/Dec25_SLEDiseaseome/Code/utils.R")

check.packages(c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr",
                 "tidyr","doParallel","caret","pbapply","BiocParallel","tibble",
                 "pathMED","NbClust","ConsensusClusterPlus","SNFtool","igraph",
                 "BloodGen3Module","purrr","vegan","SNFtool","UpSetR","ggpubr",
                 "ComplexHeatmap","tidyverse","tidytext","SnowballC","pheatmap"))

load(paste0(getwd(),"/RData/DB_step3.RData"))


## Optimize parameters based on:
#' False Discovery Rate (fdr) < 5% (Not more than 5% of significant randomly 
#' created genesets)
#' Loss of information (Variability across patients) < 50 % compared to all 
#' variability
#' Database diversity (Shannon index)
#' Additional filters: 
#' Shared genesets in at least 1/3 of total number of datasets and in at least
#' the 10% of the patients


## ······································································ Step 1 
## Get false discovery rate ----

## Create database of random pathways and get scores
all_genes <- unique(unlist(custom.db))
sizes <- c(3,quantile(sapply(custom.db, length),probs = 0.9))
random.db <- lapply(1:1000, function(x) {
  sample(all_genes, sample(sizes[1]:sizes[2], 1))})
names(random.db)<-paste("rp.",1:length(random.db))

SCORES.rd<-mScores_createReference(refObject=DATA,geneSets = random.db,cores = 12)


## Get fdr for each combination based on mscore significances
#' fdr: number of random paths that achieves significance / total nº of paths
params <- expand.grid(percentage = seq(from = 1, to = 30, by = 1), datasets = 1:15)
params$fdr<-unlist(lapply(1:nrow(params),function(p){
  drp.rp<-mScores_filterPaths(MRef=SCORES.rd,
                              min_datasets=params[p,"datasets"],
                              perc_samples=params[p,"percentage"],
                              Pcutoff=0.05,plotMetrics = FALSE)  
  res<-length(drp.rp)
  return(as.numeric((res/1000)*100))
}))
params$datasets<-as.numeric(params$datasets)

# Heatmap plot
p3<-ggplot(params, aes(x = percentage, y = datasets, fill = fdr)) +
  geom_tile(color="black") +
  labs(title = "Random-based false positives", x = "Percentage of patients",
       y = "Number of datasets",fill = "FDR") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(0, 5, 10, 50, 100)),
    limits = c(0, 100)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p3)

rm(all_genes,sizes)

save.image(paste0(getwd(),"/RData/tmpPS_1.RData"))


## ······································································ Step 2 
## Calculate the loss of information ----

#' Get nº PCAs that explained the 80% of the varianza in patients using all paths
#' 80% were selected to avoid variance explained by individual patients 
#' (i.e. paired samples or longitudinal samples from the same patients)
explainedVar<-unlist(lapply(1:length(SCORES$mscores),function(dat){
  pca <- FactoMineR::PCA(t(SCORES$mscores[[dat]]), graph = F)
  pca_eig <- as.data.frame(pca$eig)
  pca_eig <-pca_eig[pca_eig$`cumulative percentage of variance` < 80, ]
  return(nrow(pca_eig) + 1) }))
names(explainedVar)<-names(SCORES$mscores)



pcaVars<-do.call("rbind",lapply(1:nrow(params),function(i){
  #' Get nº PCAs that explained the 80% of the variance in patients using selected
  #' paths based on parameter selection (min_datasets and perc_samples)
  print(i)
  drp<-mScores_filterPaths(MRef=SCORES,
                           min_datasets=as.numeric(params[i,"datasets"]),
                           perc_samples=params[i,"percentage"],Pcutoff=0.05,
                           plotMetrics = FALSE)
  
  pca.dat<-unlist(lapply(1:length(SCORES$mscores),function(dat){
    print(dat)
    tmp<-SCORES$mscores[[dat]]
    pca <- FactoMineR::PCA(t(SCORES$mscores[[dat]][names(drp),]), graph = F)
    pca_eig <- as.data.frame(pca$eig)
    pca_eig <-pca_eig[pca_eig$`cumulative percentage of variance` < 80, ]
    return(nrow(pca_eig) + 1)
  }))
  names(pca.dat)<-names(SCORES$mscores)
  
  res<- (pca.dat / explainedVar ) * 100
  return(res)
}))
names(pcaVars)<-names(SCORES$mscores)


params$lossInfo<-as.numeric(unlist(apply(pcaVars,1,mean)))

p4<-ggplot(params, aes(x = percentage, y = datasets, fill = lossInfo)) +
  geom_tile(color="black") +
  labs(title = "Data Variance Through PCA Reduction", x = "Percentage of patients",
       y = "Number of datasets",fill = "%PCAs") +
  scale_fill_gradientn(
    colors = c("black", "white"),
    values = scales::rescale(c(0, 100)),
    limits = c(0, 100)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p4)


## Save results
save.image(paste0(getwd(),"/RData/tmpPS_1.RData"))


## ······································································ Step 3 
## Database diversity (Shannon index) ----


## Sparse all database terms
all_terms<-ann.db[ann.db$TraitID %in% names(custom.db),"term"]

all_terms <- all_terms %>% as_tibble() %>% rename(term = value) %>%
  mutate(term = str_to_lower(term)) %>% unnest_tokens(word, term) %>%
  anti_join(stop_words)  

## Words to roots
all_terms <- all_terms %>% mutate(stem = wordStem(word))  

div.res<-do.call("rbind",lapply(1:nrow(params),function(i){
  
  drp<-mScores_filterPaths(MRef=SCORES,
                           min_datasets=as.numeric(params[i,"datasets"]),
                           perc_samples=params[i,"percentage"],Pcutoff=0.05,
                           plotMetrics = FALSE)
  drp<-names(drp)
  func_terms<-ann.db[ann.db$TraitID %in% drp,"term"]
  
  clean_terms <- func_terms %>% as_tibble() %>% rename(term = value) %>%
    mutate(term = str_to_lower(term)) %>% unnest_tokens(word, term) %>%
    anti_join(stop_words)  
  
  ## Words to roots
  clean_terms <- clean_terms %>% mutate(stem = wordStem(word))  
  
  ## Here specific terms can be removed (optional)
  
  # Shannon
  freq_table <- clean_terms %>% count(stem, sort = TRUE)
  # freq_table
  term_counts <- freq_table$n
  shannon_index <- diversity(term_counts, index = "shannon")
  
  # ZCI value  
  ZCI_value <- calculate_ZCI(all_terms$stem, clean_terms$stem)
  
  res<-c("nPaths"=length(drp),"ZCI"=ZCI_value,"ShannonIndex"=shannon_index)
  return(res)
}))


params<-data.frame(params,div.res)

## Save results
save.image(paste0(getwd(),"/RData/tmpPS_1.RData"))


p5<-ggplot(params, aes(x = percentage, y = datasets, fill = ZCI)) +
  geom_tile(color="black") +
  labs(title = "ZCI", x = "Percentage of patients",
       y = "Number of datasets",fill = "ZCI") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(0.5, 1.2)),
    limits = c(0.5, 1.2)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))
plot(p5)

p5<-ggplot(params, aes(x = percentage, y = datasets, fill = ShannonIndex)) +
  geom_tile(color="black") +
  labs(title = "Shannon index", x = "Percentage of patients",
       y = "Number of datasets",fill = "Shannon index") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(3.5, 6.5)),
    limits = c(3.5, 6.5)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold")) +
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p5)



## ······································································ Step 4 
## 4: Selection of best parameters ----


## LossInfo infexion point
loss_matrix <- matrix(params$lossInfo, nrow = 30, ncol = 15)

# Get relllative change of N and P direction
grad_N <- abs(diff(loss_matrix, lag = 1, differences = 1)) 
grad_P <- abs(diff(t(loss_matrix), lag = 1, differences = 1))  


# Se relative thresholds (i.e. 10% drop)
thresholds <- quantile(params$lossInfo, probs = c(0.9, 0.8, 0.7,0.5))  

df <- params %>%
  mutate(critical_level = case_when(
    lossInfo <= thresholds[4] ~ "50% threshold",
    lossInfo <= thresholds[3] ~ "30% threshold",
    lossInfo <= thresholds[2] ~ "20% threshold",
    lossInfo <= thresholds[1] ~ "10% threshold",
    TRUE ~ "No threshold"
  ))
df$datasets<-as.numeric(df$datasets)

df_fdr<-do.call("rbind",lapply(unique(df$percentage),function(i){
  tmp<-df[df$percentage==i,]
  tmp<-tmp[tmp$fdr<=5,]
  return(tmp[1,])
}))

p6<-ggplot(df, aes(x = percentage, y = datasets, fill = ShannonIndex)) +
  geom_tile(color="black") +
  
  scale_fill_gradientn(
    colors = c("darkred","darkorange3","darkorange","gold","#86de9b"),
    values = scales::rescale(c(3.5, 6.5)),
    limits = c(3.5, 6.5)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  geom_point(data = df %>% filter(critical_level != "No threshold"), 
             aes(color = critical_level), shape = 4, size = 2) +
  scale_color_manual(values = c("10% threshold" = "#86cdeb",
                                "20% threshold" = "#619bc2",
                                "30% threshold" = "#296289",
                                "50% threshold" = "black")) +
  geom_line(data = df_fdr,aes(x = percentage, y = datasets), 
            color = "black", size = 0.3, alpha = 1) +
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1)) +
  labs(title = "Parameter optimization", x = "Percentage of patients",
       y = "Number of datasets",fill = "Shannon Index") 

plot(p6)


## Number of paths
df<-params[params$fdr<5 & params$lossInfo>50,]

p7<-ggplot(df, aes(x = percentage, y = datasets, fill = nPaths)) +
  geom_tile(color="black") +
  
  scale_fill_gradientn(
    colors = c("grey","gold","darkorange","darkorange3","darkred"),
    values = scales::rescale(c(min(df$nPaths), max(df$nPaths))),
    limits = c(min(df$nPaths), max(df$nPaths))) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1)) +
  labs(title = "Number of geneset selected", x = "Percentage of patients",
       y = "Number of datasets",fill = "DRPs") +
  geom_vline(xintercept = 10, color = "black", linetype = "dashed", size = 0.3) + 
  geom_hline(yintercept = 5, color = "black", linetype = "dashed", size = 0.3)

plot(p7)

ggarrange(p3,p4,p6,p7,nrow=2,ncol=2)


## ······································································ Step 5 
## 5: Filtering candidates ----


params.sel<-params[params$fdr<5 & params$lossInfo>50 & 
                     params$percentage>10 & params$datasets > 5,]

params.sel$names<-paste0("d",params.sel$datasets,"_",params.sel$percentage)

DRP<-lapply(1:nrow(params.sel),function(i){
  drp<-mScores_filterPaths(MRef = SCORES,
                           min_datasets=params.sel[i,"datasets"],
                           perc_samples=params.sel[i,"percentage"],
                           Pcutoff=0.05,plotMetrics = FALSE)
  # print(length(drp))
  return(names(drp))
})
names(DRP)<-params.sel$names


m1 = make_comb_mat(DRP,mode = "distinc")
UpSet(m1, set_order = params.sel$names,
      comb_order = order(comb_size(m1),decreasing = T))


rm(DRP,m1,params.sel,drp.rp,random.db,tmp,explainedVar,i,m)

save.image(paste0(getwd(),"/RData/tmpPS_1.RData"))


