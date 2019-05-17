################################################
## Results of application on BXD dataset
################################################

library(dplyr)
results_meth <- lapply(list.files("inst/extdata/BXD/",
                                  pattern="res",
                                  full.names = TRUE), readRDS)
meths <- list.files("inst/extdata/BXD/",
                      pattern="res") %>% gsub(pattern="_res.rds",replacement = "")
names(results_meth) <- meths
source("inst/R_scripts/BXD/01.loadData.R")
str(bxd)
remove_affy <- grep("Affy", colnames(bxd[[3]]))
bxd_filtered <- bxd
bxd_filtered[[3]] <- bxd_filtered[[3]][, -remove_affy]

true.clust.bxd <- rep(c("CD", "HFD"), times=c(33,31)) %>% as.factor
## Run BXD NMF
ari <- sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  tt$clust %>% mclust::adjustedRandIndex(true.clust.bxd)
})
names(ari) <- meths

xtable::xtable(ari %>% as.matrix)
## Fmeasure
fmes <- sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  tt$clust %>% FlowSOM::FMeasure(predictedClusters=true.clust.bxd %>% as.numeric() , silent = FALSE)
})
names(fmes) <- meths
xtable::xtable(cbind(fmes,ari) %>% as.matrix)


a_icluster <- results_meth[["icluster"]]$fit$beta
## Attribute names
a_icluster <- sapply(1:3, function(ii){
  rownames(a_icluster[[ii]]) <- colnames(bxd_filtered[[ii]])
  a_icluster[[ii]]
})

a_moa <- results_meth[["mocluster"]]$fit@loading
## split into a list
a_moa_tmp <- lapply(1:3, function (ll){
  a_moa[grep(paste("dat",ll, sep=""), rownames(a_moa))]
})
## Attribute names
a_moa_tmp <- sapply(1:3, function(ii){
  names(a_moa_tmp[[ii]]) <- colnames(bxd_filtered[[ii]])
  a_moa_tmp[[ii]]
})

a_sgcca <- results_meth[["sgcca"]]$fit$a

selectVars_moa <- lapply(a_moa_tmp, function (aa) aa %>% abs %>% sort(decreasing = TRUE)%>% names %>% unique  %>% head(10) )

selectVars_icluster <- lapply(a_icluster, function (aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(10)}
  )

selectVars_sgcca <- lapply(a_sgcca, function(aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(10)
  })

color <- RColorBrewer::brewer.pal(9, "Spectral")[6:8]

sapply(1:3, function (ii){
  VennDiagram::venn.diagram(list(SGCCA=selectVars_sgcca[[ii]], MoCluster=selectVars_moa[[ii]],iClusterPlus=selectVars_icluster[[ii]]) , filename=sprintf("Figs/Venn_dataBXD_%s.tiff", names(bxd)[ii]))
})
## Heatmaps



col.clust <- RColorBrewer::brewer.pal(6, "Set2")[c(4,5)]
clust_col = structure(names = c("1", "2"),col.clust)
col.true <- RColorBrewer::brewer.pal(6, "Set3")[c(5,6)]
true_col = structure(names = c("CD", "HFD"),col.true)
library(ComplexHeatmap)

ha = HeatmapAnnotation(RGCCA = results_meth[["RGCCA"]]$clust,
                       sgcca = results_meth[["sgcca"]]$clust,
                       mocluster = results_meth[["mocluster"]]$clust,
                       mcia = results_meth[["mcia"]]$clust,
                       NMF = results_meth[["NMF"]]$clust,
                       kernel = results_meth[["kernel"]]$clust,
                       icluster = results_meth[["icluster"]]$clust,
                       snf = results_meth[["snf"]]$clust,
                       col = list(RGCCA=clust_col,
                                  sgcca= clust_col,
                                  mocluster=clust_col,
                                  mcia=clust_col,
                                  NMF=clust_col,
                                  kernel=clust_col,
                                  icluster=clust_col,
                                  snf=clust_col),
                       show_legend = rep(FALSE, 9),
                       show_annotation_name = TRUE,
)
list_meth_sel <- list(icluster=selectVars_icluster,Mocluster= selectVars_moa, sgcca=selectVars_sgcca)
sapply(names(list_meth_sel), function (mm){
  print(mm)
  for( i in 1:length(bxd_filtered)){
    mat <- bxd_filtered[[i]][, list_meth_sel[[mm]][[i]]] %>% t
    f2 = circlize::colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"),
                              space = "RGB")
    if(mm=="icluster"){
      ha = HeatmapAnnotation(iClusterPlus = results_meth[[mm]]$clust,
                             Truth = true.clust.bxd,
                             col = list( iClusterPlus=clust_col, Truth=true_col),
                             show_legend = rep(FALSE, 9),
                             show_annotation_name = TRUE,
      )
    }
    if(mm=="sgcca"){
      ha = HeatmapAnnotation(SGCCA = results_meth[[mm]]$clust,
                             Truth = true.clust.bxd,
                             col = list( SGCCA=clust_col, Truth=true_col),
                             show_legend = rep(FALSE, 9),
                             show_annotation_name = TRUE,
      )
    }
    if(mm=="Mocluster"){
      ha = HeatmapAnnotation(MoCluster = results_meth[[tolower(mm)]]$clust,
                             Truth = true.clust.bxd,
                             col = list( MoCluster=clust_col, Truth=true_col),
                             show_legend = rep(FALSE, 9),
                             show_annotation_name = TRUE,
      )
    }
    ht <- Heatmap(mat,col = f2,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns=FALSE,
                  row_names_gp = gpar(fontsize =10),
                  column_names_gp = gpar(col=rep(col.true,c(33,31)), fontsize =7),
                  show_column_names = FALSE,
                  top_annotation = ha,
                  name = names(bxd_filtered)[i])
    pdf(sprintf("Figs/BXD_heatmap_%s_dataset%s.pdf", mm, names(bxd_filtered)[i]), width=8, heigh=6 )
    draw(ht,
         annotation_legend_side = "left", heatmap_legend_side = "left")
    dev.off()
  }
})

