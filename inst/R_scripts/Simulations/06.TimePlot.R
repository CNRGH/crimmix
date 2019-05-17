### Plot Time Results
library(ggplot2)
library("tibble")
library("dplyr")
pathDat <- R.utils::Arguments$getWritablePath("Time")
pathMeth <- R.utils::Arguments$getWritablePath("Time_res")
n_byClust <- sapply(1:5, function(ii) ii*c(10,20,5,25))
n <- colSums(n_byClust)
time_res <- do.call(rbind, lapply(1:4, function(ii){
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/BenchmarkTime%s", pathMeth, ii))
  print("SNF")
  ll <- list.files(pathMeth_sub, full.names = TRUE)
  rt <- do.call(rbind, lapply(ll, function(l){
    r <- readRDS(l)["elapsed",]
    data.frame(time = r, meth = gsub("_res.rds", "", gsub(paste(pathMeth_sub, "/", sep = ""),"" ,l)))
  }))
  cbind(rt, N = n[ii] %>% as.factor)
}))
time_res$meth <- gsub("NMF", "intNMF", time_res$meth)
time_res$meth <- gsub("Mocluster", "MoCluster", time_res$meth)
time_res$meth <- gsub("Kernel", "mixKernel", time_res$meth)
time_res$meth <- gsub("iCluster", "iClusterPlus", time_res$meth)


methods <- c("SNF", "RGCCA", "MCIA", "intNMF", "mixKernel","SGCCA", "MoCluster","iClusterPlus")

time_res$meth <- time_res$meth %>% factor(levels=methods)


g_time <- time_res %>% ggplot(aes(x=meth, y=log10(time)))+geom_boxplot(aes(fill=N, color=N))+scale_fill_brewer(palette="Set1")+scale_color_brewer(palette="Set1")+ylab("log(time (s))")+xlab("Methods")+theme_bw()+theme( axis.title.x = element_blank())+ guides(fill=guide_legend(title="n"), color=guide_legend(title="n"))
ggsave(g_time, filename = "Figs/Time_comp.pdf")
