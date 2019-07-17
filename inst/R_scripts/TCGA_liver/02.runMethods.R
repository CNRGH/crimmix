################################################
## Run methods on Liver cancer dataset
################################################
library(CrIMMix)
library(dplyr)
library(FlowSOM)

source("inst/R_scripts/TCGA_Liver/01.loadData.R")
str(liver)
str(clinical_data)
K <- clinical_data$pathology_T_stage %>% nlevels
### Filter
## Number of data set used for the analysis
idx_dat <- c(1,2,3,4)
set.seed(33)
## Extract patients:
path_stage <- clinical_data$pathology_T_stage %>% levels
samples <- sapply(path_stage, function(p_s){
  clinical_data$id <- rownames(clinical_data)
  sub <- clinical_data %>% filter(pathology_T_stage==p_s)
  sample(sub$id, size = 0.75*nrow(sub))
}) %>% unlist


liver_filter <- lapply(liver, function(dd){
  dd[samples,]
})
true.clust_liver <- clinical_data[samples, "pathology_T_stage"]
true.clust_liver %>% table
## filter mutations
idx.genes.non.mutated <- which(liver_filter$mutation %>% colSums==0)
liver_filter$mutation <- liver_filter$mutation[, -idx.genes.non.mutated]

SNFresults <- IntMultiOmics(liver_filter, method = "SNF", K=4)
saveRDS(SNFresults, "inst/extdata/TCGA/Liver/snf_res.rds")


c_1 <- (5/sqrt((liver_filter[idx_dat] %>% sapply(dim))[2,])) %>% round(2)
SGCCAresults <- IntMultiOmics(liver_filter[idx_dat], method = "SGCCA", K = 4, c1 = c_1,
                              ncomp = rep(5, length(liver_filter[idx_dat])))
SGCCAresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(SGCCAresults, "inst/extdata/TCGA/Liver/sgcca_res.rds")

## up to 0.139792 (bad)
MoClusterresults <- IntMultiOmics(liver_filter[idx_dat], method = "Mocluster", K = 4,
                                  k = c(0.01, 0.1, 0.01, 0.005)[idx_dat],ncomp = 4)
saveRDS(MoClusterresults, "inst/extdata/TCGA/Liver/mocluster_res.rds")
MoClusterresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)


#iClusterresults <- IntMultiOmics(liver_filter, method = "iCluster", K=K-1)
liver_filter_mcia <- liver_filter[c(1,2,4)]
MCIArresults <- IntMultiOmics(liver_filter_mcia, method = "MCIA", K=4, ncomp=5)
MCIArresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(MCIArresults, "inst/extdata/TCGA/Liver/mcia_res.rds")



MOFAresults <- IntMultiOmics(liver_filter[idx_dat], method = "MOFA", K=4)
saveRDS(MOFAresults, "inst/extdata/BXObesityD/mofa_res.rds")

Kernelresults <- IntMultiOmics(liver_filter[c(1,3)], method = "MixKernel", K=4)
Kernelresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(Kernelresults, "inst/extdata/TCGA/Liver/kernel_res.rds")

NFMresults <- IntMultiOmics(liver_filter[c(1,2)], method = "intNMF", K=4)
NFMresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(NFMresults, "inst/extdata/TCGA/Liver/NMF_res.rds")

RGCCAresults <- IntMultiOmics(liver_filter[idx_dat], method = "RGCCA", K=4, ncomp=rep(4,  length(liver_filter[idx_dat])))
RGCCAresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(RGCCAresults, "inst/extdata/TCGA/Liver/RGCCA_res.rds")

iClusterresults <- IntMultiOmics(liver_filter[idx_dat], method = "iCluster", K=3, type=c("gaussian", "binomial", "gaussian", "multinomial"))
iClusterresults$clust %>% FMeasure(predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
saveRDS(iClusterresults, "inst/extdata/TCGA/Liver/icluster_res.rds")
