library(CrIMMix)
library(dplyr)
source("inst/R_scripts/Obesity/01.loadData.R")
str(new_dat)
sapply(new_dat, dim)
true.clust_diab <- stringr::str_remove(rownames(new_dat[[1]]), "[0-9]+_")
true.clust_diab %>% table
library(mixOmics)
pca_1 <- pca(X = new_dat$meth)
pca_1 %>% plotIndiv(group = true.clust_diab)

pca_2 <- pca(X = new_dat$miRNA)
pca_2 %>% plotIndiv(group = true.clust_diab)

pca_3 <- pca(X = new_dat$mRNA, ncomp = 4)
pca_3 %>% plotIndiv(group = true.clust_diab, comp=c(1,4))


c_1 <- (2/sqrt((new_dat %>% sapply(dim))[2,])) %>% round(3)
SGCCAresults <- IntMultiOmics(new_dat, method = "SGCCA",
                              K = 2, c1 = c_1,
                              ncomp = rep(6,  length(new_dat)))
saveRDS(SGCCAresults, "inst/extdata/Obesity/sgcca_res.rds")

a <- SGCCAresults$fit$a
posDat <- lapply(a, function(aa) which(aa %>% rowSums != 0))
selectVars_sgcca <- lapply(posDat, names)
k <- c(6.556672e-05,7.279179e-04,1.037139e-03 )
MoClusterresults <- IntMultiOmics(new_dat, method="Mocluster", K=2, k= k, ncomp=6)
saveRDS(MoClusterresults, "inst/extdata/Obesity/mocluster_res.rds")

MCIArresults <- IntMultiOmics(new_dat[c(1,2,3)], method = "MCIA", K=2, ncomp=2)
saveRDS(MCIArresults, "inst/extdata/Obesity/mcia_res.rds")

SNFresults <- IntMultiOmics(new_dat[c(1,2,3)], method = "SNF", K=2, sigma=0.5, K_n=5)
saveRDS(SNFresults, "inst/extdata/Obesity/snf_res.rds")

#MOFAresults <- IntMultiOmics(new_dat, method = "MOFA", K=4)
#saveRDS(MOFAresults, "inst/extdata/BXObesityD/mofa_res.rds")

Kernelresults <- IntMultiOmics(new_dat, method = "MixKernel", K=2)
saveRDS(Kernelresults, "inst/extdata/Obesity/kernel_res.rds")

NFMresults <- IntMultiOmics(new_dat, method = "intNMF", K=2)
saveRDS(NFMresults, "inst/extdata/Obesity/NMF_res.rds")

RGCCAresults <- IntMultiOmics(new_dat, method = "RGCCA", K=2, ncomp=rep(2,  length(new_dat)))
saveRDS(RGCCAresults, "inst/extdata/Obesity/RGCCA_res.rds")

## Be careful iCluster is very slow!!! 

iClusterresults <- IntMultiOmics(new_dat, method = "iCluster", K=1,
                                 lambda=c(0.7, 0.55,0.27)*2)
iClusterresults$clust %>% table(true.clust_diab)

saveRDS(iClusterresults, "inst/extdata/Obesity/icluster_res.rds")

