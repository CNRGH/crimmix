########################
# Survival analysis 
# Performance : ARI (bad comparison with grade), but survival analysis significant
########################
library(CrIMMix)
library(CrIMMixdata)
################################################

library(dplyr)
results_meth <- lapply(list.files("inst/extdata/TCGA/liver",
                                  pattern="res",
                                  full.names = TRUE), readRDS)
meths <- list.files("inst/extdata/TCGA/liver",
                    pattern="res") %>% stringr::str_remove(pattern="_res.rds")
names(results_meth) <- meths
source("inst/R_scripts/TCGA_Liver/01.loadData.R")
str(liver)
str(clinical_data)
K <- clinical_data %>% nlevels
idx_dat <- c(1,2,3,4)
set.seed(33)
## Extract patients:
path_stage <- clinical_data$pathology_T_stage %>% levels
samples <- sapply(path_stage, function(p_s){
  clinical_data$id <- rownames(clinical_data)
  sub <- clinical_data %>% filter(pathology_T_stage==p_s)
  sample(sub$id, size = 0.75*nrow(sub))
}) %>% unlist

## Survival analysis
library(survival)
library(survminer)
`%not_in%` <- purrr::negate(`%in%`)

scrs <- list(SGCCA=do.call(cbind, results_meth[["sgcca"]]$fit$Y),
             RGCCA = do.call(cbind, results_meth[["RGCCA"]]$fit$Y),
             MCIA = results_meth[["mcia"]]$fit$mcoa$SynVar,
             SNF = results_meth[["snf"]]$fit,
             Kernel = results_meth[["kernel"]]$fit$variates$X,
             MoCluster = mogsa::moaScore(results_meth[["mocluster"]]$fit),
             icluster=results_meth[["icluster"]]$fit,
             intNMF=results_meth[["NMF"]]$fit,
             CIMLR= results_meth[["CIMLR"]]$fit,
             LRAcluster= t(results_meth[["LRAcluster"]]$fit$coordinate),
             PINSPLUS= results_meth[["PINSPLUS"]]$fit,
             ConsensusClustering= results_meth[["ConsensusClustering"]]
)
library(NbClust)
pval_5 <- sapply(names(scrs), function(ii){
  scr <- scrs[[ii]]
  print(ii)
  if(ii%in%c("icluster", "intNMF")){
    clust_surv <- scr$clusters
    names(clust_surv) = samples
  }
  else if(ii=="CIMLR"){
    clust_surv <- scr$y$cluster
    names(clust_surv) = samples
  }
  else if(ii=="PINSPLUS"){
    clust_surv <- scr$cluster1
  }
  else if(ii=="ConsensusClustering"){
    clust_surv <- scr$clust
  }
  else if(ii=="SNF"){
    K <- SNFtool::estimateNumberOfClustersGivenGraph(scr, 2:10) %>% unlist %>% table %>% which.max %>% names() %>% as.numeric
    clust_surv = scr %>% SNFtool::spectralClustering(K)
    names(clust_surv) = colnames(scr)
  }else{
    res.nbclust <- NbClust(scr, distance = "euclidean",
                           min.nc = 2, max.nc = 10, 
                           method = "ward.D2", index ="all")
    clust_surv <- res.nbclust$Best.partition
    if(sum(table(clust_surv)==1)==(length(table(clust_surv))-1)){
      t <- sapply(1:10, function(k) {
        scr %>% dist %>% hclust(method = "ward.D2") %>% cutree(k=k) %>% table
      })
      id <- which(sapply(t, function (tt) sum(tt==1)<length(tt)-1))[1]
      clust_surv <- scr %>% dist %>% hclust(method="ward.D2") %>% cutree(k=id)
    }
  }
  
  table(clust_surv)
  id_clust1 <- which(table(clust_surv)<2)
  samples_surv <- names(clust_surv)[which(clust_surv %not_in% id_clust1) ]
  data_surv <- (data.frame(clust = clust_surv[samples_surv],
                           stage = clinical_data[samples,]$pathologic_stage[samples_surv] %>% as.character,
                           time = clinical_data[samples,]$overall_survival[samples_surv],
                           status = clinical_data[samples,]$status[samples_surv]) %>% na.omit)
  
  surv_object <- Surv(time = data_surv$time %>% as.character %>% as.numeric,
                      event = data_surv$status %>% as.character %>% as.numeric)
  fit1 <- surv_fit(surv_object ~ clust, data = data_surv)
  sdf <- survdiff(surv_object ~ clust, data = data_surv)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  
  pval_res <- list(meth=ii, pval=p.val, nbByClust=paste(table(clust_surv), collapse=","))
  return(pval_res)
}) 


pval_5 



