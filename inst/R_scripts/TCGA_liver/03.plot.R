<<<<<<< HEAD
################################################
## Results of application on Liver cancer dataset
################################################

library(dplyr)
results_meth <- lapply(list.files("inst/extdata/TCGA/liver",
                                  pattern="res",
                                  full.names = TRUE), readRDS)
meths <- list.files("inst/extdata/TCGA/liver",
                      pattern="res") %>% gsub(pattern="_res.rds",replacement = "")
names(results_meth) <- meths
source("inst/R_scripts/TCGA_Liver/01.loadData.R")
str(liver)
str(clinical_data)
K <- clinical_data %>% nlevels
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
idx.genes.non.mutated <- which(liver_filter$mutation %>% colSums==0)
liver_filter$mutation <- liver_filter$mutation[, -idx.genes.non.mutated]

ari <- sapply(meths, function(m){
  print(m)
  results_meth[[m]] %>% adjustedRIComputing(true.clust_liver)
})
names(ari) <- meths

fmeas <- sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  tt$clust %>% FlowSOM::FMeasure( predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
})
names(fmeas) <- meths

xtable::xtable(cbind(fmeas, ari) %>% as.matrix)



### Explore Number of data sets used
sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  str(tt$fit)
  })
####

a_sgcca <- results_meth[["sgcca"]]$fit$a
names_dat <- lapply(a_sgcca, rownames)
posDat <- lapply(a_sgcca, function(aa) which(aa %>% rowSums != 0))


a_icluster <- results_meth[["icluster"]]$fit$beta
## Attribute names
a_icluster <- sapply(1:4, function(ii){
  rownames(a_icluster[[ii]]) <- names_dat[[ii]]
  a_icluster[[ii]]
})

a_moa <- results_meth[["mocluster"]]$fit@loading
## split into a list
a_moa_tmp <- lapply(1:4, function (ll){
  a_moa[grep(paste("dat",ll, sep=""), rownames(a_moa))]
})
## Attribute names
a_moa_tmp <- sapply(1:4, function(ii){
  names(a_moa_tmp[[ii]]) <- names_dat[[ii]]
  a_moa_tmp[[ii]]
})

selectVars_moa <- lapply(a_moa_tmp, function (aa) which(aa!=0) %>% names)

selectVars_icluster <- lapply(a_icluster, function (aa) which(rowSums(aa)!=0) %>% names)

selectVars_sgcca <- lapply(posDat, names)

cols <- RColorBrewer::brewer.pal(9,"Spectral")[c(7,8,9)]
library(grid)
### Heatmaps

a_sgcca <- results_meth[["sgcca"]]$fit$a

selectVars_moa <- lapply(a_moa_tmp, function (aa) aa %>% abs %>% sort(decreasing = TRUE)%>% names %>% unique  %>% head(20) )

selectVars_icluster <- lapply(a_icluster, function (aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(20)}
)

selectVars_sgcca <- lapply(a_sgcca, function(aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(20)
})


col.clust <- RColorBrewer::brewer.pal(6, "Set2")[1:nlevels(true.clust_liver)]
clust_col = structure(names = c("1", "2", "3","4"),col.clust)
col.true <- RColorBrewer::brewer.pal(6, "Set3")[1:nlevels(true.clust_liver)]
true_col = structure(names = levels(true.clust_liver),col.true)
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
  for( i in 1:length(liver_filter)){
    mat <- liver_filter[[i]][, list_meth_sel[[mm]][[i]]] %>% t
    f2 = circlize::colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"),
                              space = "RGB")
    ht <- Heatmap(mat,col = f2,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns=FALSE,
                  row_names_gp = gpar(fontsize =7),
                  column_names_gp = gpar(col=rep(col.true,table(true.clust_liver)), fontsize =7),
                  show_column_names = TRUE,
                  top_annotation = ha,
                  name = names(liver_filter)[i])
    pdf(sprintf("TCGA_heatmap_%s_dataset%s.pdf", mm, names(liver_filter)[i]), width=8, heigh=6 )
    draw(ht,
         annotation_legend_side = "left", heatmap_legend_side = "left")
    dev.off()
  }
})

###### Enrichissement (compute real p_values)

library(pathfindR)
RA_input <- data.frame(Gene.symbol= selectVars_moa %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_moa %>% unlist() %>% unique), 0,0.0001))
RA_output_moa <- run_pathfindR(RA_input, output_dir = "pathFindRMoCluster")
head(RA_output_moa)
RA_clustered_moa <- cluster_pathways(RA_output_moa)
RA_clustered_moa[RA_clustered$Status == "Representative", ]

RA_input_sgcca <- data.frame(Gene.symbol= selectVars_sgcca %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_sgcca %>% unlist() %>% unique), 0,0.0001))
RA_output_sgcca <- run_pathfindR(RA_input_sgcca, output_dir = "pathFindRSGCCA")
head(RA_output_sgcca)
RA_clustered <- cluster_pathways(RA_output_sgcca)
RA_clustered[RA_clustered$Status == "Representative", ]

RA_input <- data.frame(Gene.symbol= selectVars_icluster %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_icluster %>% unlist() %>% unique), 0,0.0001))
RA_output_icluster <- run_pathfindR(RA_input, output_dir = "pathFindRiCLuster")
head(RA_output_icluster)
RA_clustered <- cluster_pathways(RA_output_icluster)
RA_clustered[RA_clustered$Status == "Representative", ]




for (ii in 1:4){
vp <- VennDiagram::venn.diagram(list(Mocluster=selectVars_moa[[ii]],
                               SGCCA= selectVars_sgcca[[ii]],
                               iCluster= selectVars_icluster[[ii]]
),fill = cols,  filename = NULL, cat.cex=2, cex=2,scaled=FALSE)
pdf(sprintf("../../papers/Briefings in Bioinformatics/Figs/Venn_liver_%s.pdf", names(liver_filter)[ii]))
grid.draw(vp);
dev.off();
}

vp <- VennDiagram::venn.diagram(list(Mocluster=selectVars_moa[[2]],
                               SGCCA= selectVars_sgcca[[2]],
                               iCluster= selectVars_icluster[[2]]
), fill = cols,  filename=NULL, cat.cex=2, cex=2,scaled=FALSE)
grid.newpage()
pdf(file = "../../Poster/Fig/Venn_BXD_proteo.pdf")
grid.draw(vp);
dev.off();

vp <- VennDiagram::venn.diagram(list(Mocluster=selectVars_moa[[3]],
                               SGCCA= selectVars_sgcca[[3]],
                               iCluster= selectVars_icluster[[3]]
), fill = cols,  filename=NULL, cat.cex=2, scaled=FALSE)
grid.newpage()
  pdf(file = "../../Poster/Fig/Venn_BXD_RNA.pdf")
grid.draw(vp);
dev.off()

vp <- VennDiagram::venn.diagram(list(Mocluster=selectVars_moa[[4]],
                                     SGCCA= selectVars_sgcca[[4]],
                                     iCluster= selectVars_icluster[[4]]
), fill = cols,  filename=NULL, cat.cex=2, scaled=FALSE)
grid.newpage()
grid.draw(vp);
||||||| merged common ancestors
=======
################################################
## Results of application on Liver cancer dataset
################################################

library(dplyr)
results_meth <- lapply(list.files("inst/extdata/TCGA/liver",
                                  pattern="res",
                                  full.names = TRUE), readRDS)
meths <- list.files("inst/extdata/TCGA/liver",
                      pattern="res") %>% gsub(pattern="_res.rds",replacement = "")
names(results_meth) <- meths
source("inst/R_scripts/TCGA_Liver/01.loadData.R")
str(liver)
str(clinical_data)
K <- clinical_data %>% nlevels
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
idx.genes.non.mutated <- which(liver_filter$mutation %>% colSums==0)
liver_filter$mutation <- liver_filter$mutation[, -idx.genes.non.mutated]

ari <- sapply(meths, function(m){
  print(m)
  results_meth[[m]] %>% adjustedRIComputing(true.clust_liver)
})
names(ari) <- meths

fmeas <- sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  tt$clust %>% FlowSOM::FMeasure( predictedClusters=true.clust_liver %>% as.numeric() , silent = FALSE)
})
names(fmeas) <- meths

xtable::xtable(cbind(fmeas, ari) %>% as.matrix)



### Explore Number of data sets used
sapply(meths, function(m){
  print(m)
  tt <- results_meth[[m]]
  str(tt$fit)
  })
####

a_sgcca <- results_meth[["sgcca"]]$fit$a
names_dat <- lapply(a_sgcca, rownames)
posDat <- lapply(a_sgcca, function(aa) which(aa %>% rowSums != 0))


a_icluster <- results_meth[["icluster"]]$fit$beta
## Attribute names
a_icluster <- sapply(1:4, function(ii){
  rownames(a_icluster[[ii]]) <- names_dat[[ii]]
  a_icluster[[ii]]
})

a_moa <- results_meth[["mocluster"]]$fit@loading
## split into a list
a_moa_tmp <- lapply(1:4, function (ll){
  a_moa[grep(paste("dat",ll, sep=""), rownames(a_moa))]
})
## Attribute names
a_moa_tmp <- sapply(1:4, function(ii){
  names(a_moa_tmp[[ii]]) <- names_dat[[ii]]
  a_moa_tmp[[ii]]
})

selectVars_moa <- lapply(a_moa_tmp, function (aa) which(aa!=0) %>% names)

selectVars_icluster <- lapply(a_icluster, function (aa) which(rowSums(aa)!=0) %>% names)

selectVars_sgcca <- lapply(posDat, names)

cols <- RColorBrewer::brewer.pal(9,"Spectral")[c(7,8,9)]
library(grid)
### Heatmaps

a_sgcca <- results_meth[["sgcca"]]$fit$a

selectVars_moa <- lapply(a_moa_tmp, function (aa) aa %>% abs %>% sort(decreasing = TRUE)%>% names %>% unique  %>% head(20) )

selectVars_icluster <- lapply(a_icluster, function (aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(20)}
)

selectVars_sgcca <- lapply(a_sgcca, function(aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique %>% head(20)
})

fit_CIMLR <- results_meth[["CIMLR"]]$fit
selectVars_CIMLR <- lapply(1:4, function(ll) {
  n <- fit_CIMLR$selectfeatures$names[grep(paste("dat",ll, sep=""), fit_CIMLR$selectfeatures$names)]
  pval <- fit_CIMLR$selectfeatures$pval[grep(paste("dat",ll, sep=""), fit_CIMLR$selectfeatures$names)]
  n[which(pval<0.01)] %>% stringr::str_remove(paste("_dat",ll, sep=""))
})


col.clust <- RColorBrewer::brewer.pal(6, "Set2")[1:nlevels(true.clust_liver)]
clust_col = structure(names = c("1", "2", "3","4"),col.clust)
col.true <- RColorBrewer::brewer.pal(6, "Set3")[1:nlevels(true.clust_liver)]
true_col = structure(names = levels(true.clust_liver),col.true)
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
list_meth_sel <- list(icluster=selectVars_icluster,Mocluster= selectVars_moa, sgcca=selectVars_sgcca, CIMLR=selectVars_CIMLR)

sapply(names(list_meth_sel), function (mm){
  for( i in 1:length(liver_filter)){
    mat <- liver_filter[[i]][, list_meth_sel[[mm]][[i]]] %>% t
    f2 = circlize::colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"),
                              space = "RGB")
    ht <- Heatmap(mat,col = f2,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns=FALSE,
                  row_names_gp = gpar(fontsize =7),
                  column_names_gp = gpar(col=rep(col.true,table(true.clust_liver)), fontsize =7),
                  show_column_names = TRUE,
                  top_annotation = ha,
                  name = names(liver_filter)[i])
    pdf(sprintf("TCGA_heatmap_%s_dataset%s.pdf", mm, names(liver_filter)[i]), width=8, heigh=6 )
    draw(ht,
         annotation_legend_side = "left", heatmap_legend_side = "left")
    dev.off()
  }
})

###### Enrichissement (compute real p_values)

library(pathfindR)
RA_input <- data.frame(Gene.symbol= selectVars_moa %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_moa %>% unlist() %>% unique), 0,0.0001))
RA_output_moa <- run_pathfindR(RA_input, output_dir = "pathFindRMoCluster")
head(RA_output_moa)
RA_clustered_moa <- cluster_pathways(RA_output_moa)
RA_clustered_moa[RA_clustered$Status == "Representative", ]

RA_input_sgcca <- data.frame(Gene.symbol= selectVars_sgcca %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_sgcca %>% unlist() %>% unique), 0,0.0001))
RA_output_sgcca <- run_pathfindR(RA_input_sgcca, output_dir = "pathFindRSGCCA")
head(RA_output_sgcca)
RA_clustered <- cluster_pathways(RA_output_sgcca)
RA_clustered[RA_clustered$Status == "Representative", ]

RA_input <- data.frame(Gene.symbol= selectVars_icluster %>% unlist() %>% unique,adj.P.Val=runif(length(selectVars_icluster %>% unlist() %>% unique), 0,0.0001))
RA_output_icluster <- run_pathfindR(RA_input, output_dir = "pathFindRiCLuster")
head(RA_output_icluster)
RA_clustered <- cluster_pathways(RA_output_icluster)
RA_clustered[RA_clustered$Status == "Representative", ]




## Venn 
selectVars_moa <- lapply(a_moa_tmp, function (aa) aa %>% abs %>% sort(decreasing = TRUE)%>% names %>% unique  %>% head(100)) 

selectVars_icluster <- lapply(a_icluster, function (aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique%>% head(100)}
)

selectVars_sgcca <- lapply(a_sgcca, function(aa) {
  rowSums(aa) %>% abs %>% sort(decreasing = TRUE)  %>% names %>% unique%>% head(100)
})

fit_CIMLR <- results_meth[["CIMLR"]]$fit
selectVars_CIMLR <- lapply(1:4, function(ll) {
  n <- fit_CIMLR$selectfeatures$names[grep(paste("dat",ll, sep=""), fit_CIMLR$selectfeatures$names)]
  pval <- fit_CIMLR$selectfeatures$pval[grep(paste("dat",ll, sep=""), fit_CIMLR$selectfeatures$names)]
  n[which(pval<0.01)] %>% stringr::str_remove(paste("_dat",ll, sep=""))%>% head(100)
})


cols <- RColorBrewer::brewer.pal(10,"Spectral")[c(7,8,9,10)]

for (ii in 1:4){
vp <- VennDiagram::venn.diagram(list(Mocluster=selectVars_moa[[ii]],
                               SGCCA= selectVars_sgcca[[ii]],
                               iCluster= selectVars_icluster[[ii]],
                               CIMLR= selectVars_CIMLR[[ii]]
),fill = cols,  filename = NULL, cat.cex=2, cex=2,scaled=FALSE)
pdf(sprintf("../../papers/FigsReview/Venn_liver_%s.pdf", names(liver_filter)[ii]))
grid.draw(vp);
dev.off();
}
>>>>>>> origin/review
