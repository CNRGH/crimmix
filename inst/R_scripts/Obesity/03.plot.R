################################################
## Results of application on Obesity dataset
################################################


library(dplyr)
source("inst/R_scripts/Obesity/01.loadData.R")
results_meth <- lapply(list.files("inst/extdata/Obesity/",
                                  pattern="res",
                                  full.names = TRUE), readRDS)
meths <- list.files("inst/extdata/Obesity/",
                      pattern="res") %>% gsub(pattern="_res.rds",replacement = "")
names(results_meth) <- meths
true.clust_diab <- stringr::str_remove(rownames(new_dat[[1]]), "[0-9]+_")

ari <- sapply(meths, function(tt){
  res <- results_meth[[tt]]
  ari <- res$clust %>% mclust::adjustedRandIndex(true.clust_diab)
  return(ari)
})
names(ari) <- meths

## F-measure
fmeas <- sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  tt$clust %>% FlowSOM::FMeasure( predictedClusters=true.clust_diab %>% as.factor %>% as.numeric() , silent = FALSE)
})
names(fmeas) <- meths

xtable::xtable(cbind(fmeas, ari) %>% as.matrix)


a_icluster <- results_meth[["icluster"]]$fit$beta
## Attribute names
a_icluster <- sapply(1:3, function(ii){
  rownames(a_icluster[[ii]]) <- colnames(new_dat[[ii]])
  a_icluster[[ii]]
})

a_moa <- results_meth[["mocluster"]]$fit@loading
## split into a list
a_moa_tmp <- lapply(1:4, function (ll){
  a_moa[grep(paste("dat",ll, sep=""), rownames(a_moa))]
})
## Attribute names
a_moa_tmp <- sapply(1:3, function(ii){
  names(a_moa_tmp[[ii]]) <- colnames(new_dat[[ii]])
  a_moa_tmp[[ii]]
})

### The 10 best variables in each data set
a_sgcca <- results_meth[["sgcca"]]$fit$a

selectVars_moa <- lapply(a_moa_tmp, function (aa) aa %>% abs %>% sort(decreasing = TRUE)%>% names %>% unique  %>% head(10) )

selectVars_icluster <- lapply(a_icluster, function (aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(10)}
)

selectVars_sgcca <- lapply(a_sgcca, function(aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(10)
})

library(GeneOverlap)

for(ii in 1:3){
gom.self <- newGOM(list(iclust=selectVars_icluster[[ii]],
                        moa=selectVars_moa[[ii]],
                        sgcca=selectVars_sgcca[[ii]]),
                   genome.size=ncol(new_dat[[ii]]))
drawHeatmap(gom.self,what="Jaccard")
}


col.clust <- RColorBrewer::brewer.pal(6, "Set2")[c(4,5)]
clust_col = structure(names = c("1", "2"),col.clust)
col.true <- RColorBrewer::brewer.pal(6, "Set3")[c(5,6)]
true_col = structure(names = c("AER", "PRT"),col.true)
library(ComplexHeatmap)

id_pat <- c(grep("AER", rownames(new_dat[[1]])), grep("PRT", rownames(new_dat[[1]])))

ha = HeatmapAnnotation(RGCCA = results_meth[["RGCCA"]]$clust[id_pat],
                      sgcca = results_meth[["sgcca"]]$clust[id_pat],
                      mocluster = results_meth[["mocluster"]]$clust[id_pat],
                      mcia = results_meth[["mcia"]]$clust[id_pat],
                      NMF = results_meth[["NMF"]]$clust[id_pat],
                      kernel = results_meth[["kernel"]]$clust[id_pat],
                      icluster = results_meth[["icluster"]]$clust[id_pat],
                      snf = results_meth[["snf"]]$clust[id_pat],
                      col = list(RGCCA=clust_col,
                                 sgcca= clust_col,
                                 mocluster=clust_col,
                                 mcia=clust_col,
                                 NMF=clust_col,
                                 kernel=clust_col,
                                 icluster=clust_col,
                                 snf=clust_col
                                 ),
                      show_legend = rep(FALSE, 9),
                      show_annotation_name = TRUE,
                      )
list_meth_sel <- list(icluster=selectVars_icluster,Mocluster= selectVars_moa, sgcca=selectVars_sgcca)
sapply(names(list_meth_sel), function (mm){
  for( i in 1:length(new_dat)){
    mat <- new_dat[[i]][id_pat, list_meth_sel[[mm]][[i]]] %>% t
    f2 = circlize::colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"),
                              space = "RGB")
    ht <- Heatmap(mat,col = f2,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns=FALSE,
                  row_names_gp = gpar(fontsize =7),
                  column_names_gp = gpar(col=rep(col.true,true.clust_diab %>% table), fontsize =7),
                  show_column_names = TRUE,
                  top_annotation = ha,
                  name = names(new_dat)[i])
    pdf(sprintf("Diabete_heatmap_%s_dataset%s.pdf", mm, names(new_dat)[i]), width=8, heigh=6 )
    draw(ht,
         annotation_legend_side = "left", heatmap_legend_side = "left")
    dev.off()
  }
})


