## Influence of parameters
library(ggplot2)
library("tibble")
library("dplyr")
library(CrIMMix)
## Parallel with future
means <- c(2,2,2,2)
sds <- c(1,1,1,1)
params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means, sds, SIMPLIFY=FALSE)
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 50
nclust=4
n_byClust=c(10,20,5,25)
force=FALSE

noiseD1=c(0.2)
noiseD2=c(0.1)/10
noiseD3=c(0.1)*3
props <- c(0.005, 0.01, 0.02)

dat1 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=1000,
                  prop=props[1],params=params, noise=noiseD1)
Y1 <- dat1$data
colnames(Y1) <- sprintf("gene%s", 1:ncol(Y1))


dat2 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=500, flavor="binary",
                  params=list(c(p=0.6)), prop=props[2], noise=noiseD2)

Y2 <- dat2$data
colnames(Y2) <- sprintf("gene%s", 1:ncol(Y2))

dat3 <- simulateY(nclust=nclust,n_byClust=n_byClust, J=5000,
                  flavor="beta", params=params_beta, prop=props[3], noise=noiseD3)
Y3 <- dat3$data
colnames(Y3) <- sprintf("gene%s", 1:ncol(Y3))


sim <- list(data= list(dat1=Y1, dat2= Y2,dat3=Y3),
            biomark = list(dat1=dat1$positive,dat2=dat2$positive, dat3=dat3$positive),
            true.clust = dat1$true.clusters)


k.grid=c(0.05,0.4, 0.1)%*%t(c(0.1, 0.2,0.5,1,2,5,10))
### ncomp influence

ncomp.grid <- 1:8
Moaresults <- lapply(ncomp.grid,  IntMultiOmics, method="Mocluster",
                    data=sim$data, k=c(0.05*2, 0.05*2, 0.1), K=4 )

truth <- lapply(lapply(sim$biomark, unlist), unique)
auc_eval_moclust <- sapply (Moaresults, function(mm) {
  roc_eval(truth= truth, fit = mm$fit, method = "Mocluster")
}, simplify = FALSE)


g_moclust <- do.call(rbind, lapply(1:length(auc_eval_moclust), function (ss) {
  dd <- auc_eval_moclust[[ss]]
  n_by_data_set <- sapply(dd$TPR, length)
  tprs <- dd$TPR %>% unlist
  fprs <- dd$FPR %>% unlist
  data.frame(TPR=tprs, FPR=fprs,
             dataSet= sprintf("data set %s", rep(1:3, times=n_by_data_set)),
             k=as.factor(ss))
}))

g_moclust %>% ggplot(aes(x=FPR, y=TPR, color=k))+geom_line()+facet_grid(dataSet~.)+theme_bw()



#### SGCCA
Moaresults <- lapply(ncomp.grid,  IntMultiOmics, method="Mocluster",
                     data=sim$data, k=c(0.05*2, 0.05*2, 0.1), K=4 )

truth <- lapply(lapply(sim$biomark, unlist), unique)
auc_eval_moclust <- sapply (Moaresults, function(mm) {
  roc_eval(truth= truth, fit = mm$fit, method = "Mocluster")
}, simplify = FALSE)


g_moclust <- do.call(rbind, lapply(1:length(auc_eval_moclust), function (ss) {
  dd <- auc_eval_moclust[[ss]]
  n_by_data_set <- sapply(dd$TPR, length)
  tprs <- dd$TPR %>% unlist
  fprs <- dd$FPR %>% unlist
  data.frame(TPR=tprs, FPR=fprs,
             dataSet= sprintf("data set %s", rep(1:3, times=n_by_data_set)),
             k=as.factor(ss))
}))

g_moclust %>% ggplot(aes(x=FPR, y=TPR, color=k))+geom_line()+facet_grid(dataSet~.)+theme_bw()

