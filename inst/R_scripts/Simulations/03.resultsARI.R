## Computing ARI
source("inst/R_scripts/Simulations/00.setup.R")
pathFig <- R.utils::Arguments$getWritablePath("Figs")
library("tidyverse")
methods <- c("SNF", "RGCCA", "MCIA", "NMF", "Kernel","SGCCA", "MoCluster","iCluster")

true.clusters <- mapply(function(n, nnC) {
  rep(1:n, nnC)
}, nclust,n_by_Clust)

listBenchmark <- list.files(pathDat)


ARI_dat <- do.call(rbind, lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  
  adjustedRIs <- do.call(rbind, lapply(methods, function(mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    adjRI <- sapply(r, function(ss) {
      if (BBmisc::is.error(ss)) {
        return(NA)
      }
      if (is.na(ss)) {
        return(NA)
      } else{
        ss$clust  %>% mclust::adjustedRandIndex(true.clusters[[ii]])
      }
    })
    data.frame(method = mm, ARI = adjRI)
  }))
  adjustedRIs$noise = b
  adjustedRIs$label = sprintf("  data set 1: sd= %s,\n data set 2: prop=%s,\n  data set 3: sd=%s", noiseD1[ii], noiseD2[ii], noiseD3[ii])
  return(adjustedRIs)
}))

ARI_dat$method <- gsub("NMF", "intNMF", ARI_dat$method)
ARI_dat$method <- gsub("Kernel", "mixKernel", ARI_dat$method)



ARI_dat$method <- gsub("iCluster", "iClusterPlus", ARI_dat$method)
ARI_dat$method <- ARI_dat$method %>% factor(levels =
                                              c("SNF",
                                                "RGCCA",
                                                "MCIA",
                                                "intNMF",
                                                "mixKernel",
                                                "SGCCA",
                                                "MoCluster",
                                                "iClusterPlus"))

g_adj <- ARI_dat %>% ggplot(aes(x = method, y = ARI, fill = method)) +
  geom_boxplot() + ylab("Adjusted Rand Index") + theme_bw() +
  scale_fill_brewer(palette = "Spectral")+facet_wrap(.~noise, ncol=3, labeller = label_wrap_gen(multi_line=TRUE))+
  theme(axis.text = element_text(size = 15), strip.text.x = element_text(size = 15), legend.position = "none", axis.title.y = element_blank())+coord_flip()

g <- lemon::reposition_legend(g_adj, 'top left', panel='panel-3-3')

# define your own tags
ggsave(filename=file.path(pathFig, "Clust_eval_benchmark.pdf"),g_adj, width=13, height=13)
  
