### Integraal test
load("../data/tcga_calibvf.RDa")
data <- list(cna= cna_cal, methyl= methyl_cal, mirna=mirna_cal, rnaseq= rnaseq_cal)

library(CrIMMix)
library(dplyr)
res_Mocluster <- IntMultiOmics(data, K=2, method="Mocluster", k=c(0.1,0.8,0.6,0.9), ncomp=2)
res_snf <- IntMultiOmics(data, K=2, method="SNF", sigma=4)
true.clust <- vital_cal$vitalstatus
table(true.clust, res_Mocluster$clust)
res_Mocluster$clust %>% FlowSOM::FMeasure( predictedClusters=true.clust %>% as.numeric() , silent = FALSE)

#Recall: 0.766666666666667
res_snf$clust %>% FlowSOM::FMeasure( predictedClusters=true.clust %>% as.numeric() , silent = FALSE)
table(true.clust, res_snf$clust)


library(FactoMineR)
library("factoextra")

res.mfa <- MFA(do.call(cbind, data),
               group = sapply(data, ncol),
               name.group = names(data),
               graph = FALSE)
fviz_screeplot(res.mfa)
group <- get_mfa_var(res.mfa, "group")
group
fviz_mfa_var(res.mfa, "group")
fviz_mfa_ind(res.mfa, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
clust_AFM <- res.mfa$ind$coord %>% dist %>% hclust(method="ward.D2") %>% cutree(2)

clust_AFM %>% FlowSOM::FMeasure( predictedClusters=true.clust %>% as.numeric() , silent = FALSE)
