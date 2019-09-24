library("parallel")
source("inst/R_scripts/Simulations/00.setup.R")
listBenchmark <- list.files(pathDat)
nbCPU <- 2
for(ii in seq(along=listBenchmark)){
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
  
   print("SNF")
   SNFresults <- mclapply(data_filter, IntMultiOmics, method="SNF", K=K, mc.cores = nbCPU)
   saveRDS(SNFresults, file=file.path(pathMeth_sub, sprintf("SNF_res.rds")))
   
   print("kernel")
   Kernelresults <- mclapply(data_filter, IntMultiOmics, method="MixKernel", K=K, mc.cores = nbCPU)
   saveRDS(Kernelresults, file=file.path(pathMeth_sub, sprintf("Kernel_res.rds")))
   
   
   print("MCIA")
   MCIAresults <- mclapply(data_filter, IntMultiOmics, method="MCIA", K=K, mc.cores = nbCPU)
   saveRDS(MCIAresults, file=file.path(pathMeth_sub, sprintf("MCIA_res.rds")))
   
   
   print("Mocluster")
   Moaresults <- mclapply(data_filter, IntMultiOmics, method="Mocluster", K=K, ncomp=4, k=c(0.05,0.4, 0.1), mc.cores = nbCPU)
   saveRDS(Moaresults, file=file.path(pathMeth_sub, sprintf("Mocluster_res.rds")))
   
   print("RGCCA")
   RGCCAresults <- mclapply(data_filter, IntMultiOmics, method="RGCCA", K=K, ncomp=rep(K,3), mc.cores = nbCPU)
   saveRDS(RGCCAresults, file=file.path(pathMeth_sub, sprintf("RGCCA_res.rds")))
   
   print("NMF")
   NMFresults <- mclapply(data_filter, IntMultiOmics, method="intNMF", K=K, mc.cores = nbCPU)
   saveRDS(NMFresults, file=file.path(pathMeth_sub, sprintf("NMF_res.rds")))
   
   print("SGCCA")
   SGCCAresults <- mclapply(data_filter, IntMultiOmics, method="SGCCA", K=K, c1= c(0.3, 0.3,0.4),
                            ncomp=rep(K, 3), mc.cores = nbCPU)
   saveRDS(SGCCAresults, file=file.path(pathMeth_sub, sprintf("SGCCA_res.rds")))
   
   print("icluster")
   iCluster_results <- mclapply(data_filter, IntMultiOmics, method="iCluster", K=K-1, lambda= c(0.03, 0.03,0.03),
                                type=c("gaussian", "binomial", "gaussian"), mc.cores = nbCPU)
   saveRDS(iCluster_results, file=file.path(pathMeth_sub, sprintf("iCluster_res.rds")))
   
   
  # print("MOFA")
  # Mofaresults <- lapply(1:length(data_filter), function (dd) {
  #   IntMultiOmics(data_filter[[dd]], method="MOFA", K=K)
  # }) %>% save(file=file.path(pathMeth_sub, sprintf("MOFA_res.rds")))
  print("CIMLR")
  CIMLR_results <- lapply(data_filter, IntMultiOmics, method="CIMLR", K=K)
  saveRDS(CIMLR_results, file=file.path(pathMeth_sub, sprintf("CIMLR_res.rds")))
  
   print("LRAcluster")
   LRAcluster_results <- lapply(data_filter, IntMultiOmics, method="LRAcluster", K=K, type=c("gaussian", "binary", "gaussian"))
   saveRDS(LRAcluster_results, file=file.path(pathMeth_sub, sprintf("LRAcluster_res.rds")))
   print("PINSPLUS")
   PINSPLUS_results <- lapply(data_filter, IntMultiOmics, method="PINSPlus", K=K)
   saveRDS(PINSPLUS_results, file=file.path(pathMeth_sub, sprintf("PINSPLUS_res.rds")))
   
   print("Consensus Clustering")
   ConsensusClustering_results <- lapply(data_filter, IntMultiOmics, method="ConsensusClustering", K=K)
   saveRDS(ConsensusClustering_results, file=file.path(pathMeth_sub, sprintf("ConsensusClustering_res.rds")))
}

