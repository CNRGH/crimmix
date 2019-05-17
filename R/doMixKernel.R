#' doMixKernel
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param kernel.func the kernel function to use.
#' This parameter can be set to any user defined kernel function.
#' Widely used kernel functions are pre-implemented, that can be used by setting kernel.func
#' to one of the following strings: "kidentity", "abundance",
#' "linear" or "phylogenetic". Default: "linear".
#' @param ncomp integer. Indicates the number of components to return.
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method intNMF#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doMixKernel(data,K=4)
#' @import IntNMF
#' @import mixKernel
#' @export
doMixKernel <- function(data, K, kernel.func="linear", ncomp=2){
  dat.kernel = lapply(data, compute.kernel, kernel.func = kernel.func)
  if(is.null(names(data))){
    names(dat.kernel) <- sprintf("data_%s", 1:length(data))
  }
  if(length(data)==2){
    meta.kernel <- combine.kernels(dt1=dat.kernel[[1]], dt2=dat.kernel[[2]], method="full-UMKL")
  }
  if(length(data)==3){
    meta.kernel <- combine.kernels(dt1=dat.kernel[[1]], dt2=dat.kernel[[2]],dt3=dat.kernel[[3]], method="full-UMKL")
  }
  if(length(data)==4){
    meta.kernel <- combine.kernels(dt1=dat.kernel[[1]], dt2=dat.kernel[[2]],dt3=dat.kernel[[3]], dt4=dat.kernel[[4]],method="full-UMKL")
  }
  if(length(data)>4){
    stop("Too much blocks to run this method")
  }
  kernel.pca.result = kernel.pca(meta.kernel, ncomp = ncomp)
  clust <- kernel.pca.result$variates$X %>% dist %>% hclust (method="ward.D2") %>%
    cutree(k=K)
  res <- list(clust=clust, fit=kernel.pca.result)
}
