## ROC evaluation
source("CrIMMix/inst/Simulations_b5/R_scripts/00.setup_parameters.R")
ComputeAUC <- function(rocArray) {
  FPs_list<- sapply(1:3, function(ii) (rocArray %>% sapply(function(rr) rr$FPR))[ii,] %>% unlist %>% unique() %>% sort %>% round(2)) %>% lapply(unique)

  TPs <- do.call(rbind, lapply(c(1,2,3), function (ii){
    print(ii)
    FPs <- FPs_list[[ii]]
    tp <- sapply(FPs, FUN=function(FP) {
      sapply(1:S, FUN=function(ss) {
        roc <- rocArray[[ss]]
        ww <- which(roc$FPR[[ii]]<=FP)
        max(roc$TPR[[ii]][ww])
      })
    }) %>% colMeans
    data.frame(tp=tp, data=sprintf("dataset %s", ii))
  }))
  res <- TPs
  res$fp <- FPs_list %>% unlist
  return(res)
}
ii=5
pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, ii))

print(pathMeth_sub)
print(pathDat_sim)

list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
trueDat1 <- sapply(list.sim, function (ss) ss$biomark$dat1 %>% unlist %>% unique)
trueDat2 <-  sapply(list.sim, function (ss) ss$biomark$dat2 %>% unlist %>% unique)
trueDat3 <-  sapply(list.sim, function (ss) gsub("gene", "probe", ss$biomark$dat3 %>% unlist %>% unique) )

mm <- "Mocluster"
print(mm)
pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
ff <- readRDS(pp)

auc_eval_moclust <- sapply (1:S, function(ss) {
  f <- ff[[ss]]
  truth <- list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]])
  roc_eval(truth= truth, fit = f$fit, method = mm)
}, simplify = FALSE) %>% ComputeAUC

mm <- "SGCCA"
pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
ff <- readRDS(pp)
print(mm)

auc_eval_sgcca <- sapply (1:S, function(ss) {
  f <- ff[[ss]]
  truth <- list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]])
  roc_eval(truth= truth, fit = f$fit, method = mm)
},simplify = FALSE) %>% ComputeAUC

mm <- "iCluster"
pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
ff <- readRDS(pp)
print(mm)

auc_eval_icluster <- sapply (1:S, function(ss) {
  f <- ff[[ss]]
  truth <- list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]])
  roc_eval(truth= truth, fit = f$fit, method = mm)
},simplify = FALSE) %>% ComputeAUC

# mm <- "MOFA"
# print(mm)
#
# pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
# ff <- readRDS(pp)
# auc_eval_mofa <- sapply (1:S, function(ss) {
#   f <- ff[[ss]]
#   if(!is.na(f)){
#     truth <- list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]])
#     r <- roc_eval(truth= truth, fit = f$fit, method = mm)
#   }else{
#     r <- list(TPR = rep(0, 3),FPR = rep(0, 3) )
#   }
#   return(r)
# },simplify = FALSE) %>% ComputeAUC



ROC <- rbind(auc_eval_sgcca,
             auc_eval_moclust,
             auc_eval_icluster)
nRows <- c(nrow(auc_eval_sgcca),
           nrow(auc_eval_moclust),
           nrow(auc_eval_icluster)
)
ROC$method <- rep(c( "SGCCA", "MoCluster","icluster"),nRows)
ROC$noise = sprintf("Benchmark %s",5)

roc_eval_dat <- ROC
roc_eval_dat$method <- roc_eval_dat$method %>%
  factor(levels=c("SGCCA", "MoCluster","icluster","MOFA"))

color <- RColorBrewer::brewer.pal(9, "Spectral")[6:9]
g1 <- ggplot(roc_eval_dat, aes(x=fp, y=tp, color=method))+
  geom_line(size=1.4)+
  facet_grid(noise~data)+
  geom_abline(slope=1, intercept=0)+
  scale_color_manual(values=color)+
  theme_bw()+
  guides(fill=FALSE)+
  ylim(c(0,1))+ylab("True Positive Rate")+
  xlab("False Positive Rate")+
  theme(legend.position="bottom")
g1
ggsave(filename="Briefings in Bioinformatics/Figs/ROC_eval_benchmark_5.pdf", g1, width=7, height=4)
