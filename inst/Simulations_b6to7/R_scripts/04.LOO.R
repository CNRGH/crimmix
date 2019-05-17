pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
pathMethLOO <- R.utils::Arguments$getWritablePath("Data_Results_20181012/LOO")
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


for(ii in 1:4){
  cv <- readRDS(sprintf("%s/benchmark%s.rds", pathMethLOO,  ii))
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, ii))

  sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
  trueDat1 <- sim$biomark$dat1 %>% unlist %>% unique
  trueDat2 <-  sim$biomark$dat2 %>% unlist %>% unique
  trueDat3 <-  gsub("gene", "probe", sim$biomark$dat3 %>% unlist %>% unique)
  Truth <- list(trueDat1, trueDat2, trueDat3)
  ROC <- do.call(rbind, lapply(1:3, function (cv1) {
    mm <- methods[cv1]
    TPR <- sapply(cv[[cv1]], function (cc) sapply(roc_eval(Truth, cc$fit, method=mm)$TPR, max))
    FPR <- sapply(cv[[cv1]], function (cc) sapply(roc_eval(Truth, cc$fit, method=mm)$FPR, max))

    df <- data.frame(TPR=as.numeric(TPR),FPR=as.numeric(FPR),  meth=mm, dataset = factor(sprintf("data set %s", rep(1:3, each =60))))
  } ))

  g1 <- ROC %>% ggplot(aes(x=meth, y = TPR, colour=dataset))+geom_boxplot() + theme_bw()+ylim(c(0,1))
  g1
  ggsave(g1, filename=sprintf("Briefings in Bioinformatics/Figs/tprs_loo,benchmark%s.pdf", ii))

  g2 <- ROC %>% ggplot(aes(x=meth, y = FPR, colour=dataset))+geom_boxplot() + theme_bw()+ylim(c(0,1))
  g2
  ggsave(g2, filename=sprintf("Briefings in Bioinformatics/Figs/Fprs_loo,benchmark%s.pdf", ii))
}

# for(noise in grid.noise){
#
#   cv <- readRDS(sprintf("Data_Results_v2/LOO/noise=%s.rds", noise))
#
#   pos <- do.call(rbind, lapply(1:4, function (cv1) {
#     mm <- methods[cv1]
#     pos <- sapply(cv[[cv1]], function (cc) extractPos(fit=cc$fit, method=mm))
#     count <- do.call(rbind, lapply(c(1,2,3),function(ii){
#       pp <- pos[ii, ]
#       pp <- pp %>% unlist %>%  as.data.frame()
#       colnames(pp ) <- c("gene")
#       pp %>% group_by(gene) %>%                              # calculate the counts
#         summarize(counts = 100*n()/60) %>%
#         arrange(-counts) %>%                                # sort by counts
#         mutate(gene = factor(gene, gene), dataset=factor(sprintf("data set %s", ii)))
#     })) %>% mutate(meth=mm)
#   } ))
#
#   ### Number of genes selected at each run
#   pos %>% group_by(counts, meth, dataset) %>% summarize(n()) %>% arrange(-counts) %>% filter(counts==100) %>% xtable() %>% print()
#   ### Number of genes selected in total
#   pos %>% group_by(meth, dataset) %>% summarize(n()) %>% xtable() %>% print()
# }
#
# tpr <- function(truth, pos){
#   (pos %>% intersect(truth) %>% length)/(truth%>% length)
# }
# prec <- function (truth, pos){
#   (pos  %>% intersect(truth) %>% length)/(pos  %>% length)
# }
#
#
# F1_score_glob <- do.call(cbind, lapply(grid.noise, function(noise){
#   cv <- readRDS(sprintf("Data_Results_v2/LOO/noise=%s.rds", noise))
#   pos <- do.call(rbind, lapply(1:3, function (cv1) {
#     mm <- methods[cv1]
#     pos <- sapply(cv[[cv1]], function (cc) extractPos(fit=cc$fit, method=mm))
#     count <- do.call(rbind, lapply(c(1,2,3),function(ii){
#       pp <- pos[ii, ]
#       pp <- pp %>% unlist %>%  as.data.frame()
#       colnames(pp ) <- c("gene")
#       pp %>% group_by(gene) %>%                              # calculate the counts
#         summarize(counts = 100*n()/60) %>%
#         arrange(-counts) %>%                                # sort by counts
#         mutate(gene = factor(gene, gene), dataset=factor(sprintf("data set %s", ii)))
#     })) %>% mutate(meth=mm)
#   } )) %>% filter(counts==100)
#   pathDat_sim <- (sprintf("%s/SNR=%s", pathDat, noise))
#
#   sim <- list.files(pathDat_sim, full.names = TRUE)[[1]] %>% readRDS
#   trueDat1 <- sim$biomark$dat1 %>% unlist %>% unique
#   trueDat2 <-  sim$biomark$dat2 %>% unlist %>% unique
#   trueDat3 <-  sim$biomark$dat3 %>% unlist %>% unique
#   Truth <- list(trueDat1, trueDat2, trueDat3)
#
#   F1 <- do.call(rbind, lapply(1:3, function (ii) {
#     d <- sprintf("data set %s", ii)
#     F1_dat <- do.call(rbind,lapply(methods, function (mm){
#       pp <- pos %>% filter(meth==mm & dataset==d) %>% select(gene) %>% pull %>% as.character()
#       if(mm=="iCluster"){
#         regexp <- "[[:digit:]]+"
#         t_ii <- stringr::str_extract(Truth[[ii]], regexp) %>% as.numeric
#       }else{
#         t_ii <- Truth[[ii]]
#       }
#       F1 <- 2/(1/tpr(t_ii, pp) + 1/prec(t_ii, pp))
#     }))
#     data.frame(F1_dat)
#   }))
#   colnames(F1) <- noise
#   F1
# }))
# d <- sprintf("data set %s", 1:3)
# F1_score_glob %>% mutate(meth= rep(methods, 3), dataset=rep(d, each=4)) %>% xtable()
#
