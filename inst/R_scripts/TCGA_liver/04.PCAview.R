################################################
## PCA view Liver cancer dataset
################################################

library(CrIMMix)
library(dplyr)
source("inst/R_scripts/TCGA_Liver/01.loadData.R")
str(liver)
str(clinical_data)
K <- clinical_data$pathology_T_stage %>% nlevels
# PCA

library(ggplot2)
library(gplots)
pca_plot <- function(data) {
  data %>% ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = true.clust), show.legend = FALSE) +
    theme_bw() +
    facet_grid(. ~ title) 
}


library(mixOmics)
library(dplyr)
df.pca1 <-  pca(liver[[1]])$variates$X %>% data.frame
df.pca1$true.clust <- clinical_data$pathology_T_stage
df.pca1$title <- sprintf('Meth, Liver Cancer , PCA comp 1 - 2')
df.pca2 <-  pca(liver[[2]])$variates$X %>% data.frame
df.pca2$true.clust <-  clinical_data[,]$pathology_T_stage
df.pca2$title <- sprintf('Mutation, Liver Cancer , PCA comp 1 - 2')
df.pca3 <-  pca(liver[[3]])$variates$X %>% data.frame
df.pca3$true.clust <-  clinical_data$pathology_T_stage
df.pca3$title <- sprintf('RNA, Liver Cancer , PCA comp 1 - 2')
df.pca4 <-  pca(liver[[4]])$variates$X %>% data.frame
df.pca4$true.clust <-  clinical_data$pathology_T_stage
df.pca4$title <- sprintf('CNV, Liver Cancer , PCA comp 1 - 2')

g1 <- df.pca1 %>% pca_plot
g1
g2 <- df.pca2 %>% pca_plot
g2
g3 <- df.pca3 %>% pca_plot
g4 <- df.pca4 %>% pca_plot
g4
g <- gridExtra::grid.arrange(g1,g2,g3, g4, ncol=2)
g
