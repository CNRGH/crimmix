#' adjustedRIComputing: Compute the adjusted rand index
#'
#' @param res Output of IntMultiOmics function
#' @param true.clusters A vector containing the true clusters
#'
#' @return the Adjusted Rand Index
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data, c_2$data, c_3$data)
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="RGCCA") # ok
#' adjustedRIComputing(res,c_1$true.clusters)
#' @importFrom mclust adjustedRandIndex
#' @importFrom dplyr %>%
adjustedRIComputing <- function(res, true.clusters){
  res$clust %>% adjustedRandIndex(true.clusters)
}
