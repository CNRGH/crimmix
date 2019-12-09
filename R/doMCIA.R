<<<<<<< HEAD
#' doMCIA
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param ncomp Number of components
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method intNMF#' @export
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doMCIA(data,K=4,ncomp=3)
#' @importFrom omicade4 mcia
doMCIA <- function (data, K, ncomp=10){
  dat <- lapply(data, t)
  res <- mcia(dat, cia.nf=ncomp)
  clust <- res$mcoa$SynVar %>% dist %>% hclust(method="ward.D2") %>% cutree(k=K)
  return(list(clust=clust, fit=res))
}
||||||| merged common ancestors
=======
#' doMCIA
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param ncomp Number of components
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method MCIA
#' @export
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doMCIA(data,K=4,ncomp=3)
#' @importFrom omicade4 mcia
doMCIA <- function (data, K, ncomp=10){
  dat <- lapply(data, t)
  res <- mcia(dat, cia.nf=ncomp)
  clust <- res$mcoa$SynVar %>% dist %>% hclust(method="ward.D2") %>% cutree(k=K)
  return(list(clust=clust, fit=res))
}
>>>>>>> origin/review
