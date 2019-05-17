## PCA
library(ggplot2)
library("tibble")
library("dplyr")
library(mixOmics)
library(CrIMMix)
## Parallel with future
pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
means <- c(2,2,2,2)
sds <- c(1,1,1,1)
params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means, sds, SIMPLIFY=FALSE)
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 50
nclust=4
n_byClust=c(10,20,5,25)
force=FALSE

grid.noise <- c(0.1, 0.2, 0.5, 1)
grid.param <- list(noiseD1=c(0.2, 0.5, 0.1, 0.2),
                   noiseD2=c(0.1, 0.1, 0.5, 0.2)/10,
                   noiseD3=c(0.1, 0.1, 0.1, 0.5)*3)
props <- c(0.005, 0.01, 0.02)

for(ii in 1:4){
  nnD1 <- grid.param$noiseD1[ii]
  nnD2 <- grid.param$noiseD2[ii]
  nnD3 <- grid.param$noiseD3[ii]
  ss=1
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))
  file <- file.path(pathDat_sim, sprintf("simu%s.rds",ss ))
  dat <- readRDS(file)
  Y1 <- dat$data$dat1
  Y2 <- dat$data$dat2
  Y3 <- dat$data$dat3
  df.pca1 <-  pca(Y1)$variates$X %>% data.frame
  df.pca1$true.clust <- dat$true.clust %>% as.factor
  df.pca1$title <- sprintf('Gaussian, benchmark %s, PCA comp 1 - 2', ii)
  df.pca2 <-  pca(Y2)$variates$X %>% data.frame
  df.pca2$true.clust <- dat$true.clust %>% as.factor
  df.pca2$title <- sprintf('Binary, benchmark %s, PCA comp 1 - 2', ii)
  df.pca3 <-  pca(Y3)$variates$X %>% data.frame
  df.pca3$true.clust <- dat$true.clust %>% as.factor
  df.pca3$title <- sprintf('Beta, benchmark %s, PCA comp 1 - 2', ii)

  g1 <- df.pca1 %>% ggplot(aes(x=PC1, y=PC2))+geom_point(aes(colour=true.clust, pch=true.clust),show.legend = FALSE)+theme_bw()+facet_grid(. ~ title)+scale_color_brewer(palette="Set1")

  g2 <- df.pca2 %>% ggplot(aes(x=PC1, y=PC2))+geom_point(aes(colour=true.clust, pch=true.clust),show.legend = FALSE)+theme_bw()+facet_grid(. ~ title)+scale_color_brewer(palette="Set1")

  g3 <- df.pca3 %>% ggplot(aes(x=PC1, y=PC2))+geom_point(aes(colour=true.clust, pch=true.clust),show.legend = FALSE)+theme_bw()+facet_grid(. ~ title)+scale_color_brewer(palette="Set1")
  g <- gridExtra::grid.arrange(g1,g2,g3, ncol=3)
  ggsave(g, file=sprintf('Briefings in Bioinformatics/Figs/PCA,benchmark=%s.pdf', ii), width=15, height=7)
}
