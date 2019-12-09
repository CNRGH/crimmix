<<<<<<< HEAD
#' doiCluster
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param lambda Vector of lasso penalty terms.
#' @param type Data type, which can be gaussian, binomial, poisson, multinomial.
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method intNMF
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' \donttest{res <- doiCluster(data,K=4)}
#' @import iClusterPlus
#' @importFrom dplyr %>%
doiCluster <- function (data,  K=2, lambda = c(0.03, 0.03, 0.03, 0.03),
                         type=rep("gaussian", length(data))){
  names(data) <- sprintf("dat%s", 1:length(data))
  n_dat <- length(data)
  if(n_dat>4){
    stop('iCluster do not deal with more than 4 data sets')
  }
  if(n_dat==1){
    res <- iClusterPlus(dt1=data[[1]] ,
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==2){
    res <- iClusterPlus(dt1=data[[1]], dt2=data[[2]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==3){
    res <- iClusterPlus(dt1=data[[1]], dt2=data[[2]], dt3=data[[3]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==4){
    res <- iClusterPlus(dt1=data[[1]] , dt2=data[[2]], dt3=data[[3]], dt4=data[[4]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }

  return(list(clust=res$clusters, fit=res))
}
||||||| merged common ancestors
=======
#' doiCluster
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param lambda Vector of lasso penalty terms.
#' @param type Data type, which can be gaussian, binomial, poisson, multinomial.
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method iCluster
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' \donttest{res <- doiCluster(data,K=4)}
#' @import iClusterPlus
#' @importFrom dplyr %>%
doiCluster <- function (data,  K=2, lambda = c(0.03, 0.03, 0.03, 0.03),
                         type=rep("gaussian", length(data))){
  
  if(is.null(names(data))){
    names(data) <- sprintf("dat%s", 1:length(data))
  }
  n_dat <- length(data)
  if(n_dat>4){
    stop('iCluster do not deal with more than 4 data sets')
  }
  if(n_dat==1){
    res <- iClusterPlus(dt1=data[[1]] ,
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==2){
    res <- iClusterPlus(dt1=data[[1]], dt2=data[[2]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==3){
    res <- iClusterPlus(dt1=data[[1]], dt2=data[[2]], dt3=data[[3]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }
  if(n_dat==4){
    res <- iClusterPlus(dt1=data[[1]] , dt2=data[[2]], dt3=data[[3]], dt4=data[[4]],
                        type=type,
                        lambda=lambda,K=K,maxiter=10)
  }

  return(list(clust=res$clusters, fit=res))
}
>>>>>>> origin/review
