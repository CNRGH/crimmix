source("inst/R_scripts/Simulations/00.setup.R")
library("tidyverse")
library(lineup)
nbCPU <- 2
listBenchmark <- list.files(pathDat)
sapply( 1:8, function (ii){
  b <- listBenchmark[ii]
  K <- nclust[ii]
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))  
  ss <- 1
  sim <- list.files(pathDat_sim, full.names = TRUE)[ss] %>% readRDS
  dat <- sim$data
  comb <- list(c(1,2),c(1,3), c(2,3))
  
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
  data_filter <- dat %>% remove_zero

  t <- sapply(comb, function (cc){
    dat1 <- data_filter[cc]
    cor(dat1[[1]] , dat1[[2]] )
  })
  
})
