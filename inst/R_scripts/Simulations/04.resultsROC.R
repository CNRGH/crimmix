##########################################################################
## ROC evaluation
##########################################################################
library(dplyr)
library(CrIMMix)
library(future)
library(purrr)
plan(multiprocess)

source("inst/R_scripts/Simulations/00.setup.R")

ComputepAUC <- function(roc) {
  TPR <- roc$TPR
  FPR <- roc$FPR
  aucs <- sapply(1:length(FPR), function (ii){
    x <- FPR[[ii]]
    y <- TPR[[ii]]
    y[is.infinite(y)] <- NaN
    auc <- sum(tis::lintegrate(c(x,1), c(y,1), xint=c(x,1)))
    denom <- sum(tis::lintegrate(c(0, 0, max(x)), c(0, 1, 1), xint=c(0, 0, max(x))))
    res <- auc
    res
  })
  names(aucs) <- c("dataset 1", "dataset 2", "dataset 3")
  aucs
}

listBenchmark <- list.files(pathDat)


auc_eval_dat <- do.call(rbind,lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  
  print(pathMeth_sub)
  print(pathDat_sim)
  
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  regexp <- "[[:digit:]]+"
  
  trueDat1 <- sapply(list.sim, function (ss) ss$biomark$dat1 %>% unlist %>% unique%>% str_extract(pattern=regexp))
  trueDat2 <-  sapply(list.sim, function (ss) ss$biomark$dat2 %>% unlist %>% unique%>% str_extract(pattern=regexp))
  trueDat3 <-  sapply(list.sim, function (ss)  ss$biomark$dat3 %>% unlist %>% unique%>% str_extract(pattern=regexp))
  
  truth <- lapply(1:S, function (ss) list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]]) )
  
  
  mm <- "Mocluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  
  auc_eval_moclust <- purrr::pmap(l, roc_eval, method=mm) %>% sapply(ComputepAUC) %>% t
  
  mm <- "NMF"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_nmf <- purrr::pmap(l, roc_eval, method=mm)  %>% sapply(ComputepAUC) %>% t
  
  mm <- "MCIA"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_mcia <- purrr::pmap(l, roc_eval, method=mm) %>% sapply(ComputepAUC) %>% t
  
  mm <- "SGCCA"
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  print(mm)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_sgcca <- purrr::pmap(l, roc_eval, method=mm) %>% sapply(ComputepAUC) %>% t
  
  mm <- "RGCCA"
  pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
  ff <- readRDS(pp)
  print(mm)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_rgcca <- purrr::pmap(l, roc_eval, method=mm)  %>% sapply(ComputepAUC) %>% t
  
  
  mm <- "iCluster"
  pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
  ff <- readRDS(pp)
  print(mm)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_icluster <- purrr::pmap(l, roc_eval, method=mm)  %>% sapply(ComputepAUC) %>% t
  
  mm <- "CIMLR"
  pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
  ff <- readRDS(pp)
  print(mm)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  auc_eval_CIMLR <- purrr::pmap(l, roc_eval, method=mm)  %>% sapply(ComputepAUC) %>% t
  
  
  AUC <- rbind(auc_eval_sgcca,
               auc_eval_moclust,
               auc_eval_icluster,
               auc_eval_rgcca,
               auc_eval_nmf
               ,auc_eval_mcia, auc_eval_CIMLR)
  nRows <- c(nrow(auc_eval_sgcca),
             nrow(auc_eval_moclust),
             nrow(auc_eval_icluster),
             nrow(auc_eval_rgcca),
             nrow(auc_eval_nmf),
             nrow(auc_eval_mcia), nrow(auc_eval_CIMLR)
  )
  AUC <- AUC %>% as.data.frame %>% mutate(method = rep(c( "SGCCA", "MoCluster","icluster", "RGCCA", "iNMF", "MCIA", "CIMLR"),nRows), noise=b)
  return(AUC)
}))

saveRDS(auc_eval_dat, "inst/extdata/Data_Results_20181012/AUC_eval.rds")

auc_eval_dat <- readRDS("inst/extdata/Data_Results_20181012/AUC_eval.rds")
auc_eval_dat_2 <- auc_eval_dat %>% group_by(noise, method) %>% summarize(MeanD1= mean(`dataset 1`), MeanD2= mean(`dataset 2`), MeanD3= mean(`dataset 3`))

x1 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD1) %>% spread(method, MeanD1)
x2 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD2) %>% spread(method, MeanD2)
x3 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD3) %>% spread(method, MeanD3)
rbind(x1,x2,x3) %>% xtable::xtable()
auc_eval_dat %>% group_by(method) %>% summarize(MeanD1= mean(`dataset 1`), MeanD2= mean(`dataset 2`), MeanD3= mean(`dataset 3`)) %>% xtable::xtable()

