source("~/CrIMMix/inst/CopyOfSimulations/R_scripts/00.setup_parameters.R")
library("parallel")

pathDat_sim <- R.utils::Arguments$getReadablePath(sprintf("%s/Benchmark%s", pathDat, 5))
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, 5))

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
SNFresults <- mclapply(data, IntMultiOmics, method="SNF", K=4, mc.cores = nbCPU)
saveRDS(SNFresults, file=file.path(pathMeth_sub, sprintf("SNF_res.rds")))

print("kernel")
Kernelresults <- mclapply(data_filter, IntMultiOmics, method="MixKernel", K=4, mc.cores = nbCPU)
saveRDS(Kernelresults, file=file.path(pathMeth_sub, sprintf("Kernel_res.rds")))


print("MCIA")
MCIAresults <- mclapply(data_filter, IntMultiOmics, method="MCIA", K=4, mc.cores = nbCPU)
saveRDS(MCIAresults, file=file.path(pathMeth_sub, sprintf("MCIA_res.rds")))


print("Mocluster")
Moaresults <- mclapply(data, IntMultiOmics, method="Mocluster", K=4, ncomp=4, k=c(0.05,0.4, 0.1)*2, mc.cores = nbCPU)
saveRDS(Moaresults, file=file.path(pathMeth_sub, sprintf("Mocluster_res.rds")))

print("RGCCA")
RGCCAresults <- mclapply(data, IntMultiOmics, method="RGCCA", K=4, mc.cores = nbCPU)
saveRDS(RGCCAresults, file=file.path(pathMeth_sub, sprintf("RGCCA_res.rds")))

print("NMF")
NMFresults <- mclapply(data_filter, IntMultiOmics, method="intNMF", K=4, mc.cores = nbCPU)
saveRDS(NMFresults, file=file.path(pathMeth_sub, sprintf("NMF_res.rds")))

print("SGCCA")
SGCCAresults <- mclapply(data_filter, IntMultiOmics, method="SGCCA", K=4, c1= c(0.3, 0.3,0.4)*2,
                         ncomp=rep(3, 3), mc.cores = nbCPU)
saveRDS(SGCCAresults, file=file.path(pathMeth_sub, sprintf("SGCCA_res.rds")))

print("icluster")
iCluster_results <- mclapply(data_filter, IntMultiOmics, method="iCluster", K=3,
                             lambda= c(0.03, 0.03,0.03)*2,
                             type=c("gaussian", "binomial", "gaussian"), mc.cores = nbCPU)
saveRDS(iCluster_results, file=file.path(pathMeth_sub, sprintf("iCluster_res.rds")))

print("MOFA")
Mofaresults <- lapply(1:length(data_filter), function (dd) {
  IntMultiOmics(data_filter[[dd]], method="MOFA", K=4)
}) %>% save(file=file.path(pathMeth_sub, sprintf("MOFA_res.rds")))


