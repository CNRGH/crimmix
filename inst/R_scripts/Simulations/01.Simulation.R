### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### This script simulates benchmarks of article 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
source("inst/R_scripts/Simulations/00.setup.R")

automatic_simulation <- function(nclust, n_by_Clust,props, params_gaussian, params_beta, params_bin, noiseD1, noiseD2, noiseD3, J_gaussian, J_beta, J_binary, file, force){
  
  
  if (!file.exists(file) || force) {
    dat1 <- simulateY(nclust = nclust, n_byClust = n_by_Clust, J = J_gaussian, 
                      flavor = "normal",
                      prop = props[1], params = params_gaussian, noise = noiseD1)
    Y1 <- dat1$data
    colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))
    
    dat2 <- simulateY(nclust = nclust, n_byClust = n_by_Clust, J = J_binary, 
                      flavor = "binary",
                      params = params_bin, prop = props[2], noise = noiseD2)
    
    Y2 <- dat2$data
    colnames(Y2) <- sprintf("gene%s", 1:ncol(Y2))
    
    dat3 <- simulateY(nclust = nclust, n_byClust = n_by_Clust, J = J_beta,
                      flavor = "beta", 
                      params = list(params_beta), prop = props[3], noise = noiseD3)
    Y3 <- dat3$data
    colnames(Y3) <- sprintf("gene%s", 1:ncol(Y3))
    
    
    sim <- list(data = list(dat1 = Y1, dat2 = Y2,dat3 = Y3),
                biomark = list(dat1 = dat1$positive, dat2 = dat2$positive, dat3 = dat3$positive),
                true.clust = dat1$true.clusters)
    
    sim %>% saveRDS(file = file)
    
  }else{
    print("file already exists")
  }
  
}
force = rep(FALSE, length(nclust))
lapply(1:S, function(ss) {
  pathDat_sim <- lapply(1:length(nclust), function(bb) R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, bb)))
  file <- file.path(pathDat_sim, sprintf("simu%s.rds",ss ))
  print(file)
  l <- list(nclust, n_by_Clust,props, params_gaussian, params_beta, params_bin, noiseD1, noiseD2, noiseD3, J_gaussian, J_beta, J_binary, file, force)
  t <- purrr::pmap(l, automatic_simulation)
  
})


