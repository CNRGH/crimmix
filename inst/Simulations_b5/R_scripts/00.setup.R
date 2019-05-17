source("CrIMMix/inst/CopyOfSimulations/R_scripts/00.setup_parameters.R")
lapply (1:S, function (ss){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, 5))
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
})
