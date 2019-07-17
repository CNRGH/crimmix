#' doLRAcluster
#'
#' @param data List of matrices.
#' @param type type of data ("binary","gaussian","poisson") by default gaussian
#' @param K Number of clusters
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method LRACluster
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doSNF(data,K=4, K_n=10, sigma=0.5)
#' @import LRAcluster
#' @importFrom dplyr %>%
doLRAcluster <- function (data, K,  type= rep("gaussian", length(data))){
  dat <- lapply(data, t)
  fit <- LRAcluster(dat, dimension = K, type= type)
  clust <- fit$coordinate %>% t %>% dist %>% hclust(method="ward.D2") %>% cutree(k=K)
  res <- list(clust=clust, fit=fit)
  return(res)
}
