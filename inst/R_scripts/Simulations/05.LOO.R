##################################################################
## LOO procedure to analyse stability selection
#################################################################
library(CrIMMix)
library(dplyr)

pathMethLOO <- R.utils::Arguments$getWritablePath("inst/extdata/Data_Results_20181012_LOO")
source("inst/R_scripts/Simulations/00.setup.R")
for(ii in 1:4){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))
  
  sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
  trueDat1 <- sim$biomark$dat1 %>% unlist %>% unique
  trueDat2 <-  sim$biomark$dat2 %>% unlist %>% unique
  trueDat3 <-  sim$biomark$dat3 %>% unlist %>% unique
  Truth <- list(trueDat1, trueDat2, trueDat3)
  dat <- sim$data
  remove_zero <- function (dat){
    lapply(dat, function(dd){
      idx <- which(colSums(dd)==0)
      if(length(idx)!=0){
        return(dd[, -idx])
      }else{
        return(dd)
      }
    })
  }
  data_filter <- dat %>% remove_zero
  cv.sgcca <- LOO(data_filter, K = 4, method = "SGCCA", c1 = c(0.3, 0.3,0.4))
  cv.moclust <- LOO(data_filter, method="Mocluster", K = 4, ncomp = 4, k = c(0.05,0.4, 0.1))
  cv.icluster <- LOO(data_filter, method="iCluster", K = 3, lambda = c(0.03, 0.03,0.03))
  methods <- c("SGCCA", "Mocluster", "iCluster")
  cv <- list(cv.sgcca,cv.moclust,cv.icluster)
  saveRDS(cv, sprintf("%s/benchmark%s.rds", pathMethLOO,  ii))
}

for(ii in 1:4){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))
  print(pathDat_sim)
  sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
  trueDat1 <- sim$biomark$dat1 %>% unlist %>% unique
  trueDat2 <-  sim$biomark$dat2 %>% unlist %>% unique
  trueDat3 <-  sim$biomark$dat3 %>% unlist %>% unique
  Truth <- list(trueDat1, trueDat2, trueDat3)
  dat <- sim$data
  remove_zero <- function (dat){
    lapply(dat, function(dd){
      idx <- which(colSums(dd)==0)
      if(length(idx)!=0){
        return(dd[, -idx])
      }else{
        return(dd)
      }
    })
  }
  data_filter <- dat %>% remove_zero
  cv.CIMLR <- LOO(data_filter, K = 4, method = "CIMLR", it.max = 15)
  saveRDS(cv.CIMLR, sprintf("%s/benchmark%s_CIMLR.rds", pathMethLOO,  ii))
}

for(ii in 1:4){
  cv.CIMLR <- readRDS(sprintf("%s/benchmark%s_CIMLR.rds", pathMethLOO,  ii))
  cv <- readRDS(sprintf("%s/benchmark%s.rds", pathMethLOO,  ii))
  cv[[4]] <- cv.CIMLR
  saveRDS(cv, sprintf("%s/benchmark%s_v2.rds", pathMethLOO,  ii))
}

roc_bench <- do.call(rbind, lapply(1:4, function(b){
  cv <- readRDS(sprintf("%s/benchmark%s_v2.rds", pathMethLOO,  b))
  names(cv) <- c("SGCCA", "Mocluster", "iCluster", "CIMLR")
  counts <- do.call(rbind, lapply(c("SGCCA", "Mocluster", "iCluster", "CIMLR"), function (mm){
    cv_s <- cv[[mm]]
    print(mm)
    extract_pos_s <- sapply(cv_s, function (cv_ss){
      fit <- cv_ss$fit
      extractPos(fit, mm)
    })
    count_var <- apply(extract_pos_s, 1, function (ll){
      ll%>% unlist %>% table %>% sort(decreasing=TRUE)/length(cv_s)
    })
    
    count_var <- lapply(count_var, as.data.frame)
    n_dat <- sapply(count_var, nrow)
    count_var_dat <- do.call(rbind, count_var) %>% as.data.frame
    colnames(count_var_dat) <- c("var", "Frequency")
    count_var_dat$dataset <- sprintf("dataset %s", rep(1:3, times=n_dat))
    count_var_dat$method <- mm
    count_var_dat
  }))
  counts$benchmark <- sprintf("Benchmark %s", b)
  return(counts)
}))
library(ggplot2)
color <- RColorBrewer::brewer.pal(9, "Spectral")[6:10]
roc_bench$method <- gsub("iCluster", "iClusterPlus", roc_bench$method)
roc_bench$method <- factor(roc_bench$method, levels=c("SGCCA", "Mocluster", "iClusterPlus", "CIMLR"))
roc_bench <- roc_bench %>% mutate(data = ifelse(dataset=="dataset 1", "Gaussian",
                                         ifelse(dataset=="dataset 2", "Binary",
                                                ifelse(dataset=="dataset 3", "Beta-like", "NA")))) %>% mutate(data=factor(data, level=c("Gaussian", "Binary", "Beta-like")))


p <- ggplot(roc_bench, aes(x=method, y=Frequency, fill=method)) +
  geom_violin()+ scale_fill_manual(values=color)+facet_grid(benchmark~data)+theme_bw()+theme(legend.position = "none", axis.text = element_text(size = 15), strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15), axis.title.x = element_blank())
p
ggsave(p, filename=sprintf("../../papers/FigsReview/Loo,violin_plot_benchmark1to4.eps", b), height=10, width=15)


### Compute TPR and FPR post and pre-loo
library(xtable)
library(stringr)
listBenchmark <- sprintf("Benchmark%s", 1:8)
auc_eval_dat <- do.call(rbind,lapply(1:4, function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  
  print(pathMeth_sub)
  print(pathDat_sim)
  
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  regexp <- "[[:digit:]]+"
  
  trueDat1 <- sapply(list.sim, function (ss) ss$biomark$dat1 %>% unlist %>% unique%>% stringr::str_extract(pattern=regexp))
  trueDat2 <-  sapply(list.sim, function (ss) ss$biomark$dat2 %>% unlist %>% unique%>% stringr::str_extract(pattern=regexp))
  trueDat3 <-  sapply(list.sim, function (ss)  ss$biomark$dat3 %>% unlist %>% unique%>% stringr::str_extract(pattern=regexp))
  
  truth <- lapply(1:S, function (ss) list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]]) )
  truth1 <- truth[[1]]
  
  mm <- "Mocluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_mocluster <- ff[[1]]
  select_var_mocluster <- extractPos(all_mocluster$fit, "Mocluster")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "SGCCA"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_sgcca <- ff[[1]]
  select_var_sgcca <- extractPos(all_sgcca$fit, "SGCCA")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "iCluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_icluster <- ff[[1]]
  select_var_icluster <- extractPos(all_icluster$fit, "iCluster")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "CIMLR"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_cimlr <- ff[[1]]
  select_var_cimlr <- extractPos(all_cimlr$fit, "CIMLR")%>% lapply(str_extract, pattern=regexp)
  
  
  
  cv <- readRDS(sprintf("%s/%s_v2.rds", pathMethLOO,  tolower(b)))
  names(cv) <- c("SGCCA", "Mocluster", "iCluster", 'CIMLR')
  cv_results <- do.call(rbind, lapply(c("SGCCA", "Mocluster", "iCluster", 'CIMLR'), function(mm){
    cv_s <- cv[[mm]]
    
    extract_pos_s <- sapply(cv_s, function(cv_ss){
      fit <- cv_ss$fit
      extractPos(fit, mm)
    })
    count_var <- apply(extract_pos_s, 1, function(ll){
      ll %>% unlist %>% table %>% sort(decreasing = TRUE)/length(cv_s)
    })
    
    count_var <- lapply(count_var, as.data.frame)
    n_dat <- sapply(count_var, nrow)
    count_var_dat <- do.call(rbind, count_var) %>% as.data.frame
    colnames(count_var_dat) <- c("var", "Frequency")
    count_var_dat$dataset <- sprintf("dataset %s", rep(1:3, times = n_dat))
    count_var_dat$method <- mm
    count_var_dat
  })) %>% filter(Frequency > 0.8)
  
  stable_sgcca <- cv_results %>% filter(method =="SGCCA") %>% dplyr::select(var) %>% pull %>% as.character%>% lapply(str_extract, pattern=regexp)
  stable_mocluster <- cv_results %>% filter(method =="Mocluster") %>% dplyr::select(var) %>% pull %>% as.character%>% lapply(str_extract, pattern=regexp)
  stable_icluster <- cv_results %>% filter(method == "iCluster") %>%dplyr::select(var) %>% pull %>% as.character
  stable_cimlr <- cv_results %>% filter(method == "CIMLR") %>%dplyr::select(var) %>% pull %>% as.character%>% lapply(str_extract, pattern=regexp)
  
  sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
  
  
    fpr_func <- function(sv, Truth, sim){
    (setdiff(sv %>% unlist , Truth %>% unlist) %>% length)/(sum(sapply(sim$data, ncol)) - length(Truth %>% unlist))
  }
  tpr_func <- function(sv, Truth, sim){
    (intersect(Truth %>% unlist, sv %>% unlist) %>% length)/length(Truth %>% unlist)
  }
  
  
  fpr_before_loo <- c(sgcca = fpr_func(select_var_sgcca,truth1,sim),
                      moclust = fpr_func(select_var_mocluster,truth1,sim),
                      icluster = fpr_func(select_var_icluster,truth1,sim),
                      cimlr = fpr_func(select_var_cimlr,truth1,sim))
  
  
  fpr_after_loo <- c(sgcca = fpr_func(stable_sgcca,truth1,sim),
                     moclust = fpr_func(stable_mocluster,truth1,sim),
                     icluster = fpr_func(stable_icluster,truth1,sim),
                     cimlr = fpr_func(stable_cimlr,truth1,sim))
  
  tpr_after_loo <- c(sgcaa = tpr_func(stable_sgcca,truth1,sim),
                     moclust = tpr_func(stable_mocluster,truth1,sim),
                     icluster = tpr_func(stable_icluster,truth1,sim),
                     cimlr = tpr_func(stable_cimlr,truth1,sim)
  )
  tpr_before_loo <- c(sgcca = tpr_func(select_var_sgcca,truth1,sim),
                      moclust = tpr_func(select_var_mocluster,truth1,sim),
                      icluster = tpr_func(select_var_icluster,truth1,sim),
                      cimlr = tpr_func(select_var_cimlr,truth1,sim))
  r <- cbind(cbind(fpr_before = fpr_before_loo, fpr_after = fpr_after_loo) ,
             cbind(tpr_before = tpr_before_loo, tpr_after = tpr_after_loo)) %>% as.data.frame
  
  r$method <- rownames(r)
  r$noise <- b
  r
}))

library(tidyr)
auc_eval_dat2 <- auc_eval_dat %>%  gather(`fpr_before`, `fpr_after`, key = state_fpr, value = fpr) %>% gather(`tpr_before`, `tpr_after`, key = state_tpr, value = tpr) 
auc_eval_dat2 <- auc_eval_dat2 %>% mutate(state_fpr=factor(state_fpr, levels=c("fpr_before","fpr_after")), 
                          state_tpr=factor(state_tpr, levels=c("tpr_before","tpr_after"))) 
auc_eval_dat2$method <- gsub("icluster", "iClusterPlus", auc_eval_dat2$method)
auc_eval_dat2$method <- gsub("moclust", "MoCluster", auc_eval_dat2$method)
auc_eval_dat2$method <- gsub("sgcca", "SGCCA", auc_eval_dat2$method)
auc_eval_dat2$method <- gsub("cimlr", "CIMLR", auc_eval_dat2$method)
auc_eval_dat2$method <- factor(auc_eval_dat2$method, levels=c("SGCCA", "MoCluster", "iClusterPlus", "CIMLR"))
colnames(auc_eval_dat2)

g_fpr <- auc_eval_dat2 %>% ggplot(aes(x=method, y=fpr, fill=state_fpr))+geom_bar(stat="identity", position=position_dodge())+facet_wrap(.~noise)+theme_bw()+theme(legend.position = "top", axis.text = element_text(size = 15), strip.text.x = element_text(size = 15), axis.title.x = element_blank())+
  scale_fill_discrete(name="",
                      breaks=c("fpr_before", "fpr_after"),
                      labels=c("Raw", "Pruning by jackknife"))+ylab("False Positive Rate")
g_fpr
ggsave(g_fpr, filename = file.path("../../papers/FigsReview/FPR_LOO.pdf"), width=10, height=7)
g_tpr <- auc_eval_dat2 %>% ggplot(aes(x=method, y=tpr, fill=state_tpr))+geom_bar(stat="identity", position=position_dodge())+facet_wrap(.~noise)+theme_bw()+theme(legend.position = "top", axis.text = element_text(size = 15), strip.text.x = element_text(size = 15), axis.title.x = element_blank())+
  scale_fill_discrete(name="",
                      breaks=c("tpr_before", "tpr_after"),
                      labels=c("Raw", "Pruning by jackknife"))+ylab("True Positive Rate")

g_tpr
ggsave(g_tpr, filename = file.path("../../papers/FigsReview/TPR_LOO.pdf"), width=10, height=7)





### Count number of variables
source("inst/R_scripts/Simulations/00.setup.R")

sapply(1:4, function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))

  mm <- "Mocluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_mocluster <- ff[[1]]
  select_var_mocluster <- extractPos(all_mocluster$fit, "Mocluster")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "SGCCA"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_sgcca <- ff[[1]]
  select_var_sgcca <- extractPos(all_sgcca$fit, "SGCCA")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "iCluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_icluster <- ff[[1]]
  select_var_icluster <- extractPos(all_icluster$fit, "iCluster")%>% lapply(str_extract, pattern=regexp)
  
  mm <- "CIMLR"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  all_cimlr <- ff[[1]]
  select_var_cimlr <- extractPos(all_cimlr$fit, "CIMLR")%>% lapply(str_extract, pattern=regexp)
  
  c(mo = select_var_mocluster , 
    sgcca= select_var_sgcca ,
    iclust=select_var_icluster , 
    cimlr=select_var_cimlr )%>% sapply(length)
}) %>% xtable::xtable()
