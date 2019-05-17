#' Integration Multi-omics
#'
#' @param data List of matrices.
#' @param method A character value specifying the method to run (See Details).
#' @param ... Additionnal parameters for each methods
#'
#' @examples
#' set.seed(34)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=0.1)
#' data <- list(c_1$data, c_2$data)
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="RGCCA") # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="SGCCA") # ok
#' \dontrun{res <- CrIMMix::IntMultiOmics(data, K=4, method="MOFA") # ok}
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="SNF", K_n=10, sigma=0.5) # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="intNMF") # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="MixKernel") # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="Mocluster") # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="iCluster") # ok
#' res <- CrIMMix::IntMultiOmics(data, K=4, method="MCIA") # ok
#' @details We run the method according to the value of argument 'method':
#'
#'   If \code{method=="SGCCA"}, Variable selection for generalized
#'   canonical correlation analysis
#'   If \code{method=="RGCCA"}, Regularized generalized canonical correlation analysis.
#'   If \code{method=="iCluster"}, Integrative subtype discovery
#'   If \code{method=="intNMF"},I ntegrative clustering methods using
#'   Non-Negative Matrix Factorization
#'   If \code{method=="Mocluster"}, Identifying joint patterns across
#'   multiple omics data sets.
#'   If \code{method=="MOFA"} Multi-omics factor analysis
#'   If \code{method=="MixKernel"}, Unsupervised multiple kernel
#'   learning
#'   If \code{method=="SNF"} Run Similarity network fusion
#'   If \code{method=="MCIA"} Run Multiple Co-inertia Aalysis

#' @references
#'   Tenenhaus, A., Philippe, C., Guillemot, V., et al.
#'   Variable selection for generalized canonical correlation analysis.
#'   Biostatistics, 15(3):569–583, 2014.
#'
#'   Tenenhaus, A. and Tenenhaus, M.
#'   Regularized generalized canonical correlation analysis.
#'   Psychometrika, 76(2):257–284, 2011.
#'
#'   Shen, R., Mo, Q., Schultz, N., et al.
#'   Integrative subtype discovery in glioblastoma using icluster.
#'   PloS one, 7(4):e35236, 2012.
#'
#'   Chalise, P., Koestler, D. C., Bimali, M., et al.
#'   Integrative clustering methods for high-dimensional molecular data.
#'   Translational cancer research, 3(3):202, 2014.
#'
#'   Meng, C., Helm, D., Frejno, M., et al. mocluster:
#'   Identifying joint patterns across multiple omics data sets.
#'   Journal of proteome research, 15(3):755–765,2015.
#'
#'   Argelaguet, R., Velten, B., Arnol, D., et al.
#'   Multi-omics factor analysis - a framework for unsupervised integration
#'   of multi-omics data sets. Molecular systems biology, 14(6):e8124, 2018.
#'
#'   Mariette, J. and Villa-Vialaneix, N. Unsupervised multiple kernel
#'   learning for heterogeneous data integration. Bioinformatics, 34(6):1009–1015, 2017.
#'
#'   Wang, B., Mezlini, A. M., Demir, F., et al.
#'   Similarity network fusion for aggregating data types on a genomic scale.
#'   Nature methods, 11(3):333, 2014.
#'
#'   Meng, C., Kuster, B., Culhane, A. C., & Gholami, A. M. (2014).
#'   A multivariate approach to the integration of multi-omics datasets.
#'   BMC bioinformatics, 15(1), 162
#' @export
IntMultiOmics <- function(data, method,...){
  ## Argument
  if (!is.list(data)) {
    stop("data is not a list")
  }
  ## Check dimension
  if(max(sapply(data, dim)[1,]) != min(sapply(data, dim)[1,])){
    message(sprintf("Number of samples in dat %s is %s\n", 1:length(data), sapply(data, dim)[1,]))
    stop("data do not contains the same number of samples")
  }
  doInt <- switch(method,
                  "intNMF"=doIntNMF,
                  "iCluster"=doiCluster,
                  "MixKernel"=doMixKernel,
                  "MOFA"=doMOFA,
                  "RGCCA"= doRGCCA,
                  "SGCCA"= doSGCCA,
                  "SNF"= doSNF,
                  "Mocluster"=doMoa,
                  "MCIA"=doMCIA
  )
  res <- doInt(data,...)
  return(res)
}


