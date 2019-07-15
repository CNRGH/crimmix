#' doPINSPlus
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param sigma Variance for local model
#' @param K_n Number of nearest neighbors
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method PINSPLUS
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doSNF(data,K=4, K_n=10, sigma=0.5)
#' @import SNFtool
#' @importFrom dplyr %>%
#' @export
doPINSPlus <- function (data, K,  K_n=10, sigma=0.5){
  ## TODO !!!
  res=NA
  return(res)
}
