## Computing ARI
source("CrIMMix/inst/Simulations_b6to7/R_scripts/00.setup_parameters.R")
benchmark <- c(6,7)
ARI_dat <- do.call(rbind, lapply(1:2, function (ii){
  true.clusters <- rep(1:nclust[[ii]], n_byClust[[ii]])

  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, benchmark[[ii]]))

  adjustedRIs<- do.call(rbind, lapply(methods, function (mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    adjRI <- sapply(r, function(ss) {
      if(BBmisc::is.error(ss)){
        return(NA)
      }
      if(is.na(ss)){
        return(NA)
      } else{
        ss$clust  %>% mclust::adjustedRandIndex(true.clusters)
      }
    })
    data.frame(method=mm, ARI=adjRI)
  }))
  adjustedRIs$noise = sprintf("Benchmark %s",benchmark[[ii]])
  return(adjustedRIs)
}))
g_adj <- ARI_dat %>% ggplot(aes(x=method, y=ARI, fill=method)) +
  geom_boxplot() + ylab("Adjusted Rand Index") + theme_bw() + scale_fill_brewer(palette="Spectral") + facet_wrap(noise~., ncol=2)

g_adj
ggsave(filename=sprintf("Briefings in Bioinformatics/Figs/Clust_eval_benchmark_6_7.pdf"),g_adj, width=12, height=7)

