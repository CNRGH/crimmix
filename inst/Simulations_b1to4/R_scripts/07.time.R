## Script for time evaluation
### This script setups the parameters for the simulation
library(ggplot2)
library("tibble")
library("dplyr")
library(CrIMMix)
## Parallel with future
pathDat <- R.utils::Arguments$getWritablePath("inst/extdata/Time")
means <- c(2,2,2,2)
sds <- c(1,1,1,1)
params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means, sds, SIMPLIFY=FALSE)
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 10
nclust=4
n_byClust=sapply(1:5, function (ii) ii*c(10,20,5,25))
force=FALSE

grid.param <- list(noiseD1=c(0.2, 0.5, 0.1, 0.2),
                   noiseD2=c(0.1, 0.1, 0.5, 0.2)/10,
                   noiseD3=c(0.1, 0.1, 0.1, 0.5)*3)
props <- c(0.005, 0.01, 0.02)

nnD1 <- grid.param$noiseD1[1]
nnD2 <- grid.param$noiseD2[1]
nnD3 <- grid.param$noiseD3[1]
for(ii in 1:4){
lapply (1:S, function (ss){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/BenchmarkTime%s", pathDat, ii))
  file <- file.path(pathDat_sim, sprintf("simu%s.rds",ss ))
  print(file)
    dat1 <- simulateY(nclust=nclust,n_byClust=n_byClust[,ii], J=1000,
                      prop=props[1],params=params, noise=nnD1)
    Y1 <- dat1$data
    colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))


    dat2 <- simulateY(nclust=nclust,n_byClust=n_byClust[,ii], J=500, flavor="binary",
                      params=list(c(p=0.6)), prop=props[2], noise=nnD2)

    Y2 <- dat2$data
    colnames(Y2) <- sprintf("gene%s", 1:ncol(Y2))

    dat3 <- simulateY(nclust=nclust,n_byClust=n_byClust[,ii], J=5000,
                      flavor="beta", params=params_beta, prop=props[3], noise=nnD3)
    Y3 <- dat3$data
    colnames(Y3) <- sprintf("probe%s", 1:ncol(Y3))


    sim <- list(data= list(dat1=Y1, dat2= Y2,dat3=Y3),
                biomark = list(dat1=dat1$positive,dat2=dat2$positive, dat3=dat3$positive),
                true.clust = dat1$true.clusters)
    sim %>% saveRDS(file = file)
})
}

pathMeth <- R.utils::Arguments$getWritablePath("inst/extdata/Time_res")
for(ii in 1:4){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/BenchmarkTime%s", pathDat, ii))
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/BenchmarkTime%s", pathMeth, ii))

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
  SNFresults <- sapply(data, function(dd){
    t <- system.time(IntMultiOmics(dd, method="SNF", K=4))
  })
  saveRDS(SNFresults, file=file.path(pathMeth_sub, sprintf("SNF_res.rds")))

  print("kernel")
  Kernelresults <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd, method="MixKernel", K=4))
  })
  saveRDS(Kernelresults, file=file.path(pathMeth_sub, sprintf("Kernel_res.rds")))


  print("MCIA")
  MCIAresults <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd, method="MCIA", K=4))
  })
  saveRDS(MCIAresults, file=file.path(pathMeth_sub, sprintf("MCIA_res.rds")))


  print("Mocluster")
  Moaresults <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd,method="Mocluster", K=4, ncomp=4, k=c(0.05,0.4, 0.1)))
  })
  saveRDS(Moaresults, file=file.path(pathMeth_sub, sprintf("Mocluster_res.rds")))

  print("RGCCA")
  RGCCAresults <-sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd,method="RGCCA", K=4))
  })
  saveRDS(RGCCAresults, file=file.path(pathMeth_sub, sprintf("RGCCA_res.rds")))

  print("NMF")
  NMFresults <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd,method="intNMF", K=4))
  })
  saveRDS(NMFresults, file=file.path(pathMeth_sub, sprintf("NMF_res.rds")))

  print("SGCCA")
  SGCCAresults <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd,method="SGCCA", K=4, c1= c(0.3, 0.3,0.4),
                                   ncomp=rep(3, 3)))
  })
  saveRDS(SGCCAresults, file=file.path(pathMeth_sub, sprintf("SGCCA_res.rds")))

  print("icluster")
  iCluster_results <- sapply(data_filter, function(dd){
    t <- system.time(IntMultiOmics(dd, method="iCluster", K=3, lambda= c(0.03, 0.03,0.03),
                                   type=c("gaussian", "binomial", "gaussian")))
  })
  saveRDS(iCluster_results, file=file.path(pathMeth_sub, sprintf("iCluster_res.rds")))
}



