#' doConsensusClustering
#'
#' @param data List of matrices.
#' @param K Number of clusters
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method Consensus Clustering
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doConsensusClustering(data,K=4)
#' @import ConsensusClusterPlus
#' @importFrom dplyr %>%
doConsensusClustering <- function (data, K){
  d <- do.call(cbind, data) %>% t
  d = sweep(d,1, apply(d,1,median,na.rm=T))
  if(K==2){
    Kmax <- 10
  }else{
    Kmax <- K
  }
  fit <-  ConsensusClusterPlus(d,maxK=Kmax)
  clust <- fit[[K]]$consensusClass
  res <- list(clust=clust, fit=fit)
  return(res)
}
