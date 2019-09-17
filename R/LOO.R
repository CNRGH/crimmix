#' LOO - Integration Multi-omics
#'
#' @param data List of matrices.
#' @param method A character value specifying the method to run (See Details).
#' @param it.max An integer to define number of maximum iteration (smaller than the number of samples)
#' @param mc.cores An integer that define the number of cores to parallelize
#' @param ... Additionnal parameters for each methods
#' @return results of the method on each subset
#' @examples
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- LOO(data,K=4,method="SGCCA")
#' @export
LOO <- function(data, method, it.max=NULL,mc.cores=1, ...){
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
                  "iCluster"=doiCluster,
                  "MOFA"=doMOFA,
                  "SGCCA"= doSGCCA,
                  "Mocluster"=doMoa,
                  "CIMLR"= doCIMLR,
  )
  n <- nrow(data[[1]])
  if(is.null(it.max)){
    N <- n
  }else{
    N <- it.max
    if(N>n){
     stop(sprintf("it.max needs to be smaller than %s", n))
    }
  }

  cv <- parallel::mclapply(sample(1:N), function(i){
    cat(sprintf("fold: %s", i))
    data_i <- remove_i(data, i)
    res <- doInt(data_i, ...)
    return(res)
  }, mc.cores=mc.cores)
  return(cv)
}

