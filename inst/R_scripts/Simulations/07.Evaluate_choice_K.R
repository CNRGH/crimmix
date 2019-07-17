#### Script to evaluate the number of cluster for each method (ask by R1)

source("inst/R_scripts/Simulations/00.setup.R")
listBenchmark <- list.files(pathDat)
nbCPU <- 2
ii = 1 # Test on Benchmark 1
b <- listBenchmark[ii]
K <- nclust[ii]
pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
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
S=10
Kbest <- lapply(1:S, function(ii){
  SNF <- IntMultiOmics (data_filter[[ii]], method="SNF", K=K)
  estimationResult = SNFtool::estimateNumberOfClustersGivenGraph(SNF$fit, 2:5) %>% unlist %>% table
  K_SNF <- estimationResult %>% which.max %>% names()
  library(NbClust)
  mixKernel <- data_filter[[ii]] %>% IntMultiOmics(method="MixKernel", K=K)
  res.nbclust <- NbClust(mixKernel$fit$variates$X, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_mixKernel <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  MCIAresults <- data_filter[[ii]] %>% IntMultiOmics(method="MCIA", K=K)
  res.nbclust <- NbClust(MCIAresults$fit$mcoa$SynVar, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_MCIA <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  Moaresults <- data_filter[[ii]] %>% IntMultiOmics(method="Mocluster", K=K, ncomp=4, k=c(0.05,0.4, 0.1))
  scr <- mogsa::moaScore(Moaresults$fit)
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_Moa <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  
  
  RGCCAresults <- data_filter[[ii]] %>%IntMultiOmics( method="RGCCA", K=K)
  scr <-  do.call(cbind, RGCCAresults$fit$Y)
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_RGCCA <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  print("NMF")
  #NMFresults <- data_filter[[1]] %>%  IntMultiOmics(method="intNMF", K=K)
  ##@ TODO
  print("SGCCA")
  SGCCAresults <- data_filter[[ii]] %>% IntMultiOmics( method="SGCCA", K=K, c1= c(0.3, 0.3,0.4),
                                                      ncomp=rep(3, 3))
  scr <-  do.call(cbind, SGCCAresults$fit$Y)
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_SGCCA <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  print("icluster")
  #iCluster_results <- data_filter[[1]] %>% IntMultiOmics(method="iCluster", K=K-1, lambda= c(0.03, 0.03,0.03),
  #                             type=c("gaussian", "binomial", "gaussian"))
  ## TODO
  
  
  print("CIMLR")
  NUMC <- 2:10
  CIMLR_results <- lapply(data_filter[[ii]], t)%>%  CIMLR::CIMLR_Estimate_Number_of_Clusters(NUMC=NUMC, cores.ratio = 0)
  K_CIMLR <- NUMC[which.min(CIMLR_results$K1)]
  
  print("LRAcluster")
  LRAcluster_results <- data_filter[[ii]] %>% IntMultiOmics(method="LRAcluster", K=4, type=c("gaussian", "binary", "gaussian"))
  scr <-  LRAcluster_results$fit$coordinate %>% t
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 10, 
                         method = "ward.D2", index ="all")
  K_LRA <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  
  print("PINSPLUS")
  PINSPLUS_results <- data_filter[[ii]] %>% IntMultiOmics(method="PINSPlus", K=10)
  K_PINSPLUS <- PINSPLUS_results$clust %>% table %>% length
  
  print("Consensus Clustering")
  ConsensusClustering_results <- data_filter[[ii]] %>% IntMultiOmics(method="ConsensusClustering", K=10)
  
  K_CC <- ConsensusClustering_results$fit[[1]] %>% nrow
  K_best <- c(K_SNF, K_MCIA, K_mixKernel, K_Moa, K_SGCCA, K_RGCCA, K_CIMLR, K_LRA, K_PINSPLUS, K_CC)
  df <- cbind(K_best,method= c("SNF", "MCIA", "mixKernel", "Mocluster", "SGCCA", "RGCCA", "CIMLR", "LRACluster", "PINSPlus", "CC"))
})


df <- do.call(rbind, Kbest)%>% as.data.frame()
df <- df %>% mutate(K_best= K_best %>% as.integer) 
df <- df %>% mutate(method=as.character(method))
df$method[df$method=="CC"] <- "ConsensusClustering"
df <- df %>% mutate(method= factor(method, levels =c("SNF",
                                              "RGCCA",
                                              "MCIA",
                                              "intNMF",
                                              "mixKernel",
                                              "SGCCA",
                                              "Mocluster",
                                              "iClusterPlus",
                                              "CIMLR",
                                              "LRACluster",
                                              "PINSPlus",
                                              "ConsensusClustering")))
library(ggplot2)
 g_kbest <- df %>% ggplot(aes(x=method, y=K_best)) + geom_boxplot()+theme_bw()+geom_hline(yintercept=4, col="red", lty="dashed")+theme(legend.position = "none")+coord_flip()
 g_kbest
ggsave("../../papers/FigsReview/K_best.eps", g_kbest, width=10, height=5)
