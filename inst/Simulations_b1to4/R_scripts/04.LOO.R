pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
pathMethLOO <- R.utils::Arguments$getWritablePath("Data_Results_20181012/LOO")
pathMeth <- R.utils::Arguments$getWritablePath("Data_Results_20181012")

S <- 50
nclust=4
n_byClust=c(10,20,5,25)
nbCPU=15

grid.noise <- c(0.1, 0.2, 0.5, 1)
grid.param <- list(noiseD1=c(0.2, 0.5, 0.1, 0.2),
                   noiseD2=c(0.1, 0.1, 0.5, 0.2)/10,
                   noiseD3=c(0.1, 0.1, 0.1, 0.5)*3)
props <- c(0.005, 0.01, 0.02)

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

library(CrIMMix)
library(dplyr)
roc_bench <- do.call(rbind, lapply(1:4, function(b){
  cv <- readRDS(sprintf("%s/benchmark%s.rds", pathMethLOO,  b))
  names(cv) <- c("SGCCA", "Mocluster", "iCluster")
  test <- do.call(rbind, lapply(c("SGCCA", "Mocluster", "iCluster"), function (mm){
    cv_s <- cv[[mm]]

    extract_pos_s<- sapply(cv_s, function (cv_ss){
      fit <- cv_ss$fit
      extractPos(fit, mm)
    })
    count_var <- apply(extract_pos_s, 1, function (ll){
      ll%>% unlist %>% table %>% sort(decreasing=TRUE)/60 %>% as.numeric
    })

    count_var <- lapply(count_var, as.data.frame)
    n_dat <- sapply(count_var, nrow)
    count_var_dat <- do.call(rbind, count_var) %>% as.data.frame
    colnames(count_var_dat) <- c("var", "Frequency")
    count_var_dat$dataset <- sprintf("dataset %s", rep(1:3, times=n_dat))
    count_var_dat$method <- mm
    count_var_dat
  }))
  test$benchmark <- sprintf("Benchmark %s", b)
  return(test)
}))
library(ggplot2)
color <- RColorBrewer::brewer.pal(9, "Spectral")[6:9]

roc_bench$method <- factor(roc_bench$method, levels=c("SGCCA", "Mocluster", "iCluster"))
p <- ggplot(roc_bench, aes(x=method, y=Frequency, fill=method)) +
  geom_violin()+ scale_fill_manual(values=color)+facet_grid(benchmark~dataset)+theme_bw()
p
ggsave(p, filename=sprintf("../../papers/Briefings in Bioinformatics/Figs/Loo,violin_plot_benchmark1to4.pdf", b), height=10, width=15)


### Compute TPR and FPR post and pre-loo
library(xtable)
res_post_loo <- lapply(1:4, function(b){
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, b))
  sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
  trueDat1 <- sim$biomark$dat1 %>% unlist %>% unique
  trueDat2 <-  sim$biomark$dat2 %>% unlist %>% unique
  trueDat3 <-  gsub("gene", "probe", sim$biomark$dat3 %>% unlist %>% unique)
  Truth <- list(trueDat1, trueDat2, trueDat3)
  cv <- readRDS(sprintf("%s/benchmark%s.rds", pathMethLOO,  b))
  names(cv) <- c("SGCCA", "Mocluster", "iCluster")
  test <- do.call(rbind, lapply(c("SGCCA", "Mocluster", "iCluster"), function (mm){
    cv_s <- cv[[mm]]

    extract_pos_s<- sapply(cv_s, function (cv_ss){
      fit <- cv_ss$fit
      extractPos(fit, mm)
    })
    count_var <- apply(extract_pos_s, 1, function (ll){
      ll%>% unlist %>% table %>% sort(decreasing=TRUE)/60 %>% as.numeric
    })

    count_var <- lapply(count_var, as.data.frame)
    n_dat <- sapply(count_var, nrow)
    count_var_dat <- do.call(rbind, count_var) %>% as.data.frame
    colnames(count_var_dat) <- c("var", "Frequency")
    count_var_dat$dataset <- sprintf("dataset %s", rep(1:3, times=n_dat))
    count_var_dat$method <- mm
    count_var_dat
  }))
  pathRes <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathMeth, b))
  all_sgcca <- readRDS(file.path(pathRes, "SGCCA_res.rds"))[[1]]
  select_var_sgcca <- extractPos(all_sgcca$fit, "SGCCA")
  all_mocluster <- readRDS(file.path(pathRes, "Mocluster_res.rds"))[[1]]
  select_var_mocluster<- extractPos(all_mocluster$fit, "Mocluster")
  all_iCluster <- readRDS(file.path(pathRes, "iCluster_res.rds"))[[1]]
  select_var_icluster <- lapply(extractPos(all_iCluster$fit, "iCluster"), function(l) paste("gene",l, sep=""))

  (setdiff(select_var_sgcca %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist))
  (setdiff(select_var_mocluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist))
  (setdiff(select_var_icluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist))


  stable_sgcca <- test %>% filter(Frequency>0.8 & method=="SGCCA") %>% select(var) %>% pull %>% as.character
  stable_mocluster <- test %>% filter(Frequency>0.8 & method=="Mocluster") %>% select(var) %>% pull %>% as.character
  stable_icluster <- test %>% filter(Frequency>0.8 & method=="iCluster") %>% mutate(var=paste("gene",var, sep="")) %>% select(var) %>% pull %>% as.character

  fpr_before_loo <- c(sgcca=(setdiff(select_var_sgcca %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)),
  moclust= (setdiff(select_var_mocluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)),
  icluster= (setdiff(select_var_icluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)))


  fpr_after_loo <-c(sgcca=(setdiff(stable_sgcca %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)),
  moclust=(setdiff(stable_mocluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)),
  icluster=(setdiff(stable_icluster %>% unlist , Truth %>% unlist)%>% length)/(sum(sapply(sim$data, ncol))-length(Truth%>% unlist)))

  tpr_after_loo <- c(sgcaa=(intersect(Truth %>% unlist, stable_sgcca) %>% length)/length(Truth %>% unlist),
  mclust=(intersect(Truth %>% unlist, stable_mocluster) %>% length)/length(Truth %>% unlist),
  icluster=(intersect(Truth %>% unlist, stable_icluster) %>% length)/length(Truth %>% unlist)
)
  tpr_before_loo <- c(sgcca=(intersect(Truth %>% unlist, select_var_sgcca %>% unlist) %>% length)/length(Truth%>% unlist),
  moclust=(intersect(Truth %>% unlist, select_var_mocluster %>% unlist) %>% length)/length(Truth %>% unlist),
  icluster=(intersect(Truth %>% unlist, select_var_icluster %>% unlist) %>% length)/length(Truth %>% unlist))
  print("FPR")
  r <- cbind(cbind(fpr_before= fpr_before_loo, fpr_after= fpr_after_loo) ,
    cbind(tpr_before= tpr_before_loo, tpr_after= tpr_after_loo)) %>% as.data.frame
})


