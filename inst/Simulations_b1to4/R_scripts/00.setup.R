### This script setups the parameters for the simulation
library(ggplot2)
library("tibble")
library("dplyr")
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
  lapply (1:S, function (ss){
    pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))
    file <- file.path(pathDat_sim, sprintf("simu%s.rds",ss ))
    print(file)
    if(!file.exists(file)||force){
      dat1 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=1000,
                        prop=props[1],params=params, noise=nnD1)
      Y1 <- dat1$data
      colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))


      dat2 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=500, flavor="binary",
                        params=list(c(p=0.6)), prop=props[2], noise=nnD2)

      Y2 <- dat2$data
      colnames(Y2) <- sprintf("gene%s", 1:ncol(Y2))

      dat3 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=5000,
                        flavor="beta", params=params_beta, prop=props[3], noise=nnD3)
      Y3 <- dat3$data
      colnames(Y3) <- sprintf("probe%s", 1:ncol(Y3))


      sim <- list(data= list(dat1=Y1, dat2= Y2,dat3=Y3),
                  biomark = list(dat1=dat1$positive,dat2=dat2$positive, dat3=dat3$positive),
                  true.clust = dat1$true.clusters)
      sim %>% saveRDS(file = file)
    }else{
      print("file already exists")
    }
  }
  )
}
