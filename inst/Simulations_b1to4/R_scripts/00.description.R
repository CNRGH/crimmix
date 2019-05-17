## Additionnal analysis
library(ggplot2)
library("tibble")
library("dplyr")
library(CrIMMix)
library(mixOmics)
## Parallel with future
means <- c(2,2,2,2)
sds <- c(1,1,1,1)
params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means, sds, SIMPLIFY=FALSE)
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 50
nclust=4
n_byClust=c(10,20,5,25)
force=FALSE

grid.noise <- c(0.1, 0.2, 0.5, 1)

for(noise in grid.noise){
  dat1 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=1000,
                    prop=0.005,params=params, noise=noise*2)
  Y1 <- dat1$data
  colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))
  pdf(sprintf('Briefings in Bioinformatics/Figs/PCA-Gaussian,noise=%s.pdf', noise))
  pca(Y1) %>% plotIndiv(group=dat1$true.clusters %>% as.factor,
                        title = sprintf('Gaussian, noise=%s, PCA comp 1 - 2', noise))
  dev.off()
}

for(noise in grid.noise){
  dat2 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=5000,
                    flavor="beta", params=params_beta, prop=0.02, noise=noise*3)
  Y2 <- dat2$data
  colnames(Y2) <- sprintf("probe%s", 1:ncol(Y2))
  pdf(sprintf('Briefings in Bioinformatics/Figs/PCA-Beta,noise=%s.pdf', noise))
  pca(Y2) %>% plotIndiv(group=dat2$true.clusters %>% as.factor,
                        title = sprintf('Beta, noise=%s, PCA comp 1 - 2', noise))
  dev.off()
}

for(noise in grid.noise){
  dat3 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=500, flavor="binary",
                    params=list(c(p=0.6)), prop=0.01, noise=noise/10)

  Y3 <- dat3$data
  colnames(Y3) <- sprintf("gene%s", 1:ncol(Y3))
  pdf(sprintf('Briefings in Bioinformatics/Figs/PCA-Binary,noise=%s.pdf', noise))
  pca(Y3) %>% plotIndiv(group=dat3$true.clusters %>% as.factor,
                        title = sprintf('Binary, noise=%s, PCA comp 1 - 2', noise))
  dev.off()
}
