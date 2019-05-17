library(dplyr)
library(CrIMMix)
## ROC evaluation
pathDat <- R.utils::Arguments$getReadablePath("Data_sim_20181012")
pathMeth <- R.utils::Arguments$getReadablePath("Data_Results_20181012")
S <- 50
nclust=4
n_byClust=c(10,20,5,25)
nbCPU=15

grid.noise <- c(0.1, 0.2, 0.5, 1)
grid.param <- list(noiseD1=c(0.2, 0.5, 0.1, 0.2),
                   noiseD2=c(0.1, 0.1, 0.5, 0.2)/10,
                   noiseD3=c(0.1, 0.1, 0.1, 0.5)*3)
props <- c(0.005, 0.01, 0.02)

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

roc_eval_dat <- do.call(rbind, lapply(1:4, function(ii){
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

  mm <- "NMF"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)

  auc_eval_nmf <- sapply (1:S, function(ss) {
    f <- ff[[ss]]
    truth <- list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]])
    roc_eval(truth= truth, fit = f$fit, method = mm)
  }, simplify = FALSE) %>% ComputeAUC

  mm <- "MCIA"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)

  auc_eval_mcia <- sapply (1:S, function(ss) {
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

  mm <- "RGCCA"
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  print(mm)

  auc_eval_rgcca <- sapply (1:S, function(ss) {
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
               auc_eval_icluster,
               auc_eval_rgcca,
               auc_eval_nmf
               ,auc_eval_mcia)
  nRows <- c(nrow(auc_eval_sgcca),
             nrow(auc_eval_moclust),
             nrow(auc_eval_icluster),
             nrow(auc_eval_rgcca),
             nrow(auc_eval_nmf),
             nrow(auc_eval_mcia)
  )
  ROC$method <- rep(c( "SGCCA", "MoCluster","icluster", "RGCCA", "iNMF", "MCIA"),nRows)
  ROC$noise = sprintf("Benchmark %s",ii)
  return(ROC)
}))

library(ggplot2)
roc_eval_dat %>% head
roc_eval_dat$method <- roc_eval_dat$method %>%
  factor(levels=c("SNF", "RGCCA", "MCIA", "iNMF", "Kernel","SGCCA", "MoCluster","icluster"))
roc_eval_dat$relevantBM <- NA
roc_eval_dat$relevantBM[which(roc_eval_dat$data=='dataset 1')] <- "Relevant biomarkers = 0.5%"
roc_eval_dat$relevantBM[which(roc_eval_dat$data=='dataset 2')] <- "Relevant biomarkers = 1%"
roc_eval_dat$relevantBM[which(roc_eval_dat$data=='dataset 3')] <- "Relevant biomarkers = 2%"

color <- RColorBrewer::brewer.pal(8, "Spectral")[c(2,3,4,6,7,8)]
g1 <- ggplot(roc_eval_dat, aes(x=fp, y=tp, color=method))+
  geom_line(size=0.7)+
  facet_grid(noise~data+relevantBM, )+
  geom_abline(slope=1, intercept=0)+
  scale_color_manual(values=color)+
  theme_bw()+
  guides(fill=FALSE)+
  ylim(c(0,1))+ylab("True Positive Rate")+
  xlab("False Positive Rate")+
  theme(legend.position="top", strip.text = element_text(face="bold", size=11))
g1

ggsave(filename="../../Poster/Fig/ROC_eval_benchmark.pdf", g1, width=7, height=7)
