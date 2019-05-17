## Computing ARI
pathFig <- R.utils::Arguments$getWritablePath("Figs")
library("tidyverse")
methods <- c("SNF", "RGCCA", "MCIA", "NMF", "Kernel","SGCCA", "MoCluster","iCluster")
pathMeth <- R.utils::Arguments$getReadablePath("Data_Results_20181012")
nclust=4
n_byClust=c(10,20,5,25)
force=FALSE

grid.noise <- c(0.1, 0.2, 0.5, 1)
grid.param <- list(noiseD1=c(0.2, 0.5, 0.1, 0.2),
                   noiseD2=c(0.1, 0.1, 0.5, 0.2)/10,
                   noiseD3=c(0.1, 0.1, 0.1, 0.5)*3)
props <- c(0.005, 0.01, 0.02)

true.clusters <- rep(1:nclust, n_byClust)

ARI_dat <- do.call(rbind, lapply(1:4, function (ii){
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, ii))

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
  adjustedRIs$noise = sprintf("Benchmark %s",ii)
  return(adjustedRIs)
}))

ARI_dat$method <- gsub("NMF", "iNMF", ARI_dat$method)
ARI_dat$method <- ARI_dat$method %>% factor(levels=
                                              c("SNF",
                                                "RGCCA",
                                                "MCIA",
                                                "iNMF",
                                                "Kernel",
                                                "SGCCA",
                                                "MoCluster",
                                                "iCluster"))
g_adj <- ARI_dat %>% ggplot(aes(x=method, y=ARI, fill=method)) +
  geom_boxplot() + ylab("Adjusted Rand Index") + theme_bw() +
  scale_fill_brewer(palette="Spectral") +facet_wrap(. ~ noise, ncol=2) +
  theme(legend.position="none",axis.text=element_text(size=12))

dat_text <- data.frame(
  label = c("data set 1: sd=0.2\n data set 2: prop=0.01\n data set 3: sd=0.3",
            "data set 1: sd=0.5\n data set 2: prop=0.01\n data set 3: sd=0.3",
            "data set 1: sd=0.1\n data set 2: prop=0.05\n data set 3: sd=0.3",
            "data set 1: sd=0.2\n data set 2: prop=0.02\n data set 3: sd=1.5"),
  noise = c( "Benchmark 1",  "Benchmark 2",  "Benchmark 3",  "Benchmark 4"),
  method=c("SNF",
           "RGCCA",
           "MCIA",
           "iNMF",
           "Kernel",
           "SGCCA",
           "MoCluster",
           "iCluster"))

g_adj2 <- g_adj+ geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -0.2,
  size=5
)+theme(strip.text = element_text(face="bold", size=16), axis.title = element_text(size=20))
g_adj2
ggsave(filename=sprintf("../../Poster/Fig/Clust_eval_benchmark.pdf"),g_adj2, width=13, height=8)

