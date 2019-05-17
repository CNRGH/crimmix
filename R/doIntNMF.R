#' doIntNMF
#'
#' @param data List of matrices.
#' @param K Number of latent variables and number of clusters
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method intNMF
#' @export
#'
#' @examples
#' set.seed(33)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' \donttest{res <- doIntNMF(data,K=4)}
#' @import IntNMF
#' @importFrom dplyr %>%
doIntNMF <- function (data, K){
  dat <- lapply(data, function (dd){
    if (!all(dd>=0)) dd <- pmax(dd + abs(min(dd)), .Machine$double.eps)
    dd <- dd/max(dd)
    return(dd %>% as.matrix)
  })

  # The function nmf.mnnals requires the samples to be on rows and variables on columns.
  result.intNMF <- dat %>% IntNMF::nmf.mnnals(k=K)
  clust.intNMF <- result.intNMF$clusters
  res <- list(clust=clust.intNMF, fit = result.intNMF)
  return(res)
}
