################################################
## Run methods of application on BXD dataset
################################################
library(CrIMMix)
library(dplyr)
source("inst/R_scripts/BXD/01.loadData.R")
str(bxd)
remove_affy <- grep("Affy", colnames(bxd[[3]]))
bxd_filtered <- bxd
bxd_filtered[[3]] <- bxd_filtered[[3]][, -remove_affy]
true.clust.bxd <- rep(c("CD", "HFD"), each=32)
c_1 <- (10/sqrt((bxd_filtered %>% sapply(dim))[2,])) %>% round(2)
SGCCAresults <- IntMultiOmics(bxd_filtered, method = "SGCCA", K=2, c1=c_1)
saveRDS(SGCCAresults, "inst/extdata/BXD/sgcca_res.rds")
a <- SGCCAresults$fit$a
posDat <- lapply(a, function(aa) which(aa %>% rowSums != 0))
selectVars_sgcca <- lapply(posDat, names)

k <- selectVars_sgcca %>% sapply(length)/(bxd_filtered %>% sapply(dim))[2,]
MoClusterresults <- IntMultiOmics(bxd_filtered, method = "Mocluster", K=2, k=k)
saveRDS(MoClusterresults, "inst/extdata/BXD/mocluster_res.rds")



library(iClusterPlus)
## Run on bioinfo (slow slow slow!!!)
# cv.fit = iClusterPlus::tune.iClusterPlus(cpus=10,
#                                          dt1=bxd_filtered[[1]],
#                                          dt2=bxd_filtered[[2]],
#                            dt3=bxd_filtered[[3]],
#                            type=c("gaussian","gaussian","gaussian"),K=1,n.lambda=35,
#                            scale.lambda=c(1,1,1),maxiter=20)
cv.fit <- readRDS("inst/extdata/BXD/iclustCV.rds")
ii <- sapply(cv.fit$fit, function (dd) dd$clusters %>% mclust::adjustedRandIndex(true.clust.bxd)) %>% which.max
iClusterresults <- IntMultiOmics(bxd_filtered, method = "iCluster", K=1,
                                 lambda=cv.fit$lambda[ii,])
saveRDS(iClusterresults, "inst/extdata/BXD/icluster_res.rds")

MCIArresults <- IntMultiOmics(bxd_filtered, method = "MCIA", K=2)
saveRDS(MCIArresults, "inst/extdata/BXD/mcia_res.rds")

SNFresults <- IntMultiOmics(bxd_filtered, method = "SNF", K=2)
saveRDS(SNFresults, "inst/extdata/BXD/snf_res.rds")

MOFAresults <- IntMultiOmics(bxd_filtered, method = "MOFA", K=2)
saveRDS(MOFAresults, "inst/extdata/BXD/mofa_res.rds")

Kernelresults <- IntMultiOmics(bxd_filtered[c(1,3)], method = "MixKernel", K=2)
saveRDS(Kernelresults, "inst/extdata/BXD/kernel_res.rds")

NFMresults <- IntMultiOmics(bxd_filtered[c(1,3)], method = "intNMF", K=2)
saveRDS(NFMresults, "inst/extdata/BXD/NMF_res.rds")

RGCCAresults <- IntMultiOmics(bxd_filtered, method = "RGCCA", K=2)
saveRDS(RGCCAresults, "inst/extdata/BXD/RGCCA_res.rds")


