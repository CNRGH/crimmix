source("CrIMMix/inst/Simulations_b6to7/R_scripts/00.setup_parameters.R")
benchmark <- c(6,7)
for(ii in 1:2){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, benchmark[ii]))
  lapply (1:S, function (ss){
    params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means[[ii]], sds[[ii]], SIMPLIFY=FALSE)

    file <- file.path(pathDat_sim, sprintf("simu%s.rds",ss ))
    print(file)
    if(!file.exists(file)||force){
      dat1 <- simulateY(nclust=nclust[[ii]],n_byClust=n_byClust[[ii]], J=1000,
                        prop=props[1],params=params, noise=nnD1)
      Y1 <- dat1$data
      colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))


      dat2 <- simulateY(nclust=nclust[[ii]],n_byClust=n_byClust[[ii]], J=500, flavor="binary",
                        params=list(c(p=0.6)), prop=props[2], noise=nnD2)

      Y2 <- dat2$data
      colnames(Y2) <- sprintf("gene%s", 1:ncol(Y2))

      dat3 <- simulateY(nclust=nclust[[ii]],n_byClust=n_byClust[[ii]], J=5000,
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

  })

}
