## Combine Datasets
source("inst/R_scripts/Simulations/00.setup.R")
library("tidyverse")
methods <- c("SNF","CIMLR", "LRAcluster", "PINSPLUS", "ConsensusClustering", "RGCCA", "MCIA", "NMF", "Kernel","SGCCA", "MoCluster","iCluster")
listBenchmark <- list.files(pathDat)
nbCPU <- 2
ii <- 1
b <- listBenchmark[ii]
K <- nclust[ii]
pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
pathMeth_comb <- "inst/extdata/Data_results_20190822_comb"
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth_comb, b))
print(pathMeth_sub)
list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
data <- lapply(list.sim, function (ll) ll$data)
remove_zero <- function (dat){
  lapply(dat, function(dd){
    idx <- which(colSums(dd)==0)
    if(length(idx)!=0){
      return(dd[, -idx])
    }else{
      return(dd)
    }
  })
}
data_filter <- data %>% lapply(remove_zero)
data1 <- data_filter[[1]]
print("SNF")

comb <- list(c(1:3), 1:2, 2:3, c(1,3))
lapply(1:10, function (rr){
  data1 <- data_filter[[rr]]
  res <- lapply(1:length(comb), function (cc){
    c <- comb[[cc]]
    SNFresults <- IntMultiOmics(data1[c], method="SNF", K=K)
    #saveRDS(SNFresults, file=file.path(pathMeth_sub, sprintf("SNF_res_comb%s.rds", cc)))
    
    print("kernel")
    Kernelresults <- IntMultiOmics(data1[c],  method="MixKernel", K=K)
    #saveRDS(Kernelresults, file=file.path(pathMeth_sub, sprintf("Kernel_res_comb%s.rds", cc)))
    
    print("MCIA")
    MCIAresults <- IntMultiOmics(data1[c],  method="MCIA", K=K)
    #saveRDS(MCIAresults, file=file.path(pathMeth_sub, sprintf("MCIA_res_comb%s.rds", cc)))
    
    print("Mocluster")
    Moaresults <-  IntMultiOmics(data1[c],  method="Mocluster", K=K, ncomp=4, k=c(0.05,0.4, 0.1)[c])
    #saveRDS(Moaresults, file=file.path(pathMeth_sub, sprintf("Mocluster_res_comb%s.rds", cc)))
    
    print("RGCCA")
    RGCCAresults <-  IntMultiOmics(data1[c],  method="RGCCA", K=K)
    #saveRDS(RGCCAresults, file=file.path(pathMeth_sub, sprintf("RGCCA_res_comb%s.rds", cc)))
    
    print("NMF")
    NMFresults <- IntMultiOmics(data1[c],  method="intNMF", K=K)
    #saveRDS(NMFresults, file=file.path(pathMeth_sub, sprintf("NMF_res_comb%s.rds", cc)))
    
    print("SGCCA")
    SGCCAresults <-  IntMultiOmics(data1[c],  method="SGCCA", K=K, c1= c(0.3, 0.3,0.4)[c],
                                   ncomp=rep(3, length(c)))
    #saveRDS(SGCCAresults, file=file.path(pathMeth_sub, sprintf("SGCCA_res_comb%s.rds", cc)))
    
    print("icluster")
    iCluster_results <-  IntMultiOmics(data1[c],  method="iCluster", K=K-1, lambda= c(0.03, 0.03,0.03)[c],
                                       type=c("gaussian", "binomial", "gaussian")[c])
    #saveRDS(iCluster_results, file=file.path(pathMeth_sub, sprintf("iCluster_res_comb%s.rds", cc)))
    
    print("CIMLR")
    CIMLR_results <-  IntMultiOmics(data1[c],  method="CIMLR", K=K)
    #saveRDS(CIMLR_results, file=file.path(pathMeth_sub, sprintf("CIMLR_res_comb%s.rds", cc)))
    
    print("LRAcluster")
    LRAcluster_results <- IntMultiOmics(data1[c],  method="LRAcluster", K=K, type=c("gaussian", "binary", "gaussian")[c])
    #saveRDS(LRAcluster_results, file=file.path(pathMeth_sub, sprintf("LRAcluster_res_comb%s.rds", cc)))
    print("PINSPLUS")
    PINSPLUS_results <-  IntMultiOmics(data1[c],  method="PINSPlus", K=K)
    #saveRDS(PINSPLUS_results, file=file.path(pathMeth_sub, sprintf("PINSPLUS_res_comb%s.rds", cc)))
    #
    print("Consensus Clustering")
    ConsensusClustering_results <-  IntMultiOmics(data1[c],  method="ConsensusClustering", K=K)
    #saveRDS(ConsensusClustering_results, file=file.path(pathMeth_sub, sprintf("ConsensusClustering_res_comb%s.rds", cc)))
    res_tot <- list(SNFresults$clust, Kernelresults$clust, MCIAresults$clust, Moaresults$clust, NMFresults$clust, RGCCAresults$clust, SGCCAresults$clust, iCluster_results$clust, CIMLR_results$clust, LRAcluster_results$clust, PINSPLUS_results$clust, ConsensusClustering_results$clust)
    names(res_tot) <- c("SNF","mixKernel","MCIA","MoCluster",
                        "intNMF","RGCCA","SGCCA","iClusterPlus","CIMLR",
                        "LRAcluster","PINSPLUS","ConsensusClustering")
    return(res_tot)
  })
  names(res) <- sprintf("com%s", 1:4)
  saveRDS(res, file=file.path(pathMeth_sub, sprintf("res_sim%s.rds", rr)))
})





true.clusters <- mapply(function(n, nnC) {
  rep(1:n, nnC)
}, nclust,n_by_Clust)


ARI_dat <- do.call(rbind, lapply(1:10, function(ii){
  files <- list.files(pathMeth_sub, pattern=sprintf("sim%s.rds", ii), full.names = TRUE)
  fil1 <- readRDS(files)
  dat <- do.call(cbind, lapply(fil1, function (r){
    adjRI <- r  %>% sapply(mclust::adjustedRandIndex,true.clusters[[1]])
  }))
  dat <- dat %>% data.frame %>% mutate(method=rownames(dat), sim=sprintf("sim%s", ii)) %>% gather(com1:com4, value=ARI, key=comb)
  return(dat)
}))

comb_list <- comb
ARI_dat <- ARI_dat %>% mutate(comb = as.numeric(as.factor(comb)))

combination=sapply(ARI_dat$comb, function (cc) paste(comb_list[[cc]], collapse=","))

methods <- c("SNF",
             "LRAcluster",
             "PINSPLUS",
             "ConsensusClustering",
             "RGCCA",
             "MCIA",
             "intNMF",
             "mixKernel",
             "SGCCA",
             "MoCluster",
             "iClusterPlus",
             "CIMLR"
)

ARI_dat <- ARI_dat %>% mutate(combination=factor(combination, levels=sapply(comb_list, paste,collapse=",")))


g_ari <- ARI_dat %>% ggplot(aes(y=ARI, x=combination, fill=factor(method, levels=methods)))+  geom_boxplot() + facet_wrap(.~method, ncol=3)+scale_fill_brewer(palette = "Set3")+theme_bw()+theme(legend.position = "none")


ggsave(filename=file.path("../../papers/FigsReview/", "combination.eps"),g_ari, width=13, height=13)


