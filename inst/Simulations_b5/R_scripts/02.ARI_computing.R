## Computing ARI
source("CrIMMix/inst/Simulations_b5/R_scripts/00.setup_parameters.R")
methods <- c("SNF", "RGCCA", "MCIA", "NMF", "Kernel","SGCCA", "MoCluster","iCluster")
true.clusters <- rep(1:nclust, n_byClust)
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, 5))

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
adjustedRIs$noise <- "Benchmark 5"
g_adj <- adjustedRIs %>% ggplot(aes(x=method, y=ARI, fill=method)) +
  geom_boxplot() + ylab("Adjusted Rand Index") + theme_bw() + scale_fill_brewer(palette="Spectral") + facet_wrap(noise~., ncol=1)

g_adj
ggsave(filename=sprintf("Briefings in Bioinformatics/Figs/Clust_eval_benchmark_5.pdf"),g_adj, width=12, height=7)

