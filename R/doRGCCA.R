#' doRGCCA
#'
#' @param data List of matrices.
#' @param K Nmber of clusters
#' @param C A design matrix that describes the relationships between blocks (default: complete design).
#' @param ncomp A 1 * J vector that contains the numbers of c
#' omponents for each block (default: \code{rep(1, length(data}),
#' which gives one component per block.)
#' @param scheme The value is "horst", "factorial", "centroid"
#' or any diffentiable convex scheme function g designed by
#' the user (default: "centroid").
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
#' res <- doRGCCA(data,K=4)
#' @importFrom RGCCA rgcca
#' @export
doRGCCA <- function (data, K, C=1-diag(length(data)),
                     ncomp=rep(1, length(data)), scheme="centroid"){
  tau = rep(1,length(data))
  ## rgcca algorithm using the dual formulation for X1 and X2
  ## and the dual formulation for X3
  result.rgcca = data %>% rgcca(C, tau, ncomp = ncomp,
                                scheme = scheme,
                                verbose = FALSE)
  resDat <- do.call(cbind, result.rgcca$Y)
  clust.rgcca <- resDat %>% dist %>% hclust(method="ward.D2") %>% cutree(k=K)
  res <- list(clust=clust.rgcca, fit = result.rgcca)
  return(res)
}
