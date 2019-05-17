#' doMoa
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param ncomp An integer; the number of components to calculate.
#' To calculate more components requires longer computational time.
#' @param k The absolute number (if k >= 1) or the proportion (if 0<k<1) of non-zero coefficients
#'  for the variable loading vectors. It could be a single value or
#'  a vector has the same length as x so the sparsity of individual matrix could be different.
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
#' res <- doMoa(data,K=4)
#' @import mogsa
#' @importFrom dplyr %>%
#' @export
doMoa <- function (data, K, ncomp=2, k=0.05){
  names(data) <- sprintf("dat%s", 1:length(data))
  data_t <- lapply(data, t)
  moas <- data_t %>% mbpca(ncomp = ncomp, k =k, method = "globalScore", option = "lambda1",
                        center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000,verbose=FALSE)
  scrs <- moas %>% moaScore
  clust.moa <- scrs %>% dist %>% hclust(method="ward.D2") %>% cutree(K)
  res <- list(clust=clust.moa, fit=moas)
}
