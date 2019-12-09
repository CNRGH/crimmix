<<<<<<< HEAD
#' doMOFA
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method MOFA
#'
#' @examples
#' \dontrun{c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doMOFA(data,K=4)}
#' @importFrom dplyr %>%
#' @import MOFAtools
doMOFA <- function(data, K){
  data_t <- lapply(data, t)
  MOFAobject <- createMOFAobject(data_t)
  TrainOptions <- getDefaultTrainOptions()
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  ModelOptions$numFactors <- 20
  DataOptions <- getDefaultDataOptions()
  TrainOptions$DropFactorThreshold <- 0.1
  n_inits <- 4
  MOFAlist <- lapply(1:n_inits, function(it) {

    TrainOptions$seed <- 2018 + it

    MOFAobject <- prepareMOFA(
      MOFAobject,
      DataOptions = DataOptions,
      ModelOptions = ModelOptions,
      TrainOptions = TrainOptions
    )
    out <- tryCatch(
      {
        runMOFA(MOFAobject)
      },
      error=function(cond) {
        message(paste("Shut down all components, no structure found in the data."))
        # Choose a return value in case of error
        return(NA)
      }
    )
    return(out)
  })
  ## all is na
  if(sum(sapply(MOFAlist, is.na))==4){
    res <- NA
  }else{
    if(sum(sapply(MOFAlist, is.na))>=1){
      MOFAlist <- MOFAlist[-which(is.na(MOFAlist))]
    }
    MOFAobject <- selectModel(MOFAlist, plotit = FALSE)
    clust.MOFA <- MOFAobject@Expectations$Z%>%
      dist %>% hclust(method="ward.D2") %>% cutree(k=K)
    res <- list(clust=clust.MOFA, fit = MOFAobject)
  }
  return(res)

}
||||||| merged common ancestors
=======
#' doMOFA
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method MOFA
#'
#' @examples
#' \dontrun{c_1 <- simulateY(J=100, prop=0.1, noise=1)
#' c_2 <- simulateY(J=200, prop=0.1, noise=1)
#' c_3 <- simulateY(J=50, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doMOFA(data,K=4)}
#' @importFrom dplyr %>%
doMOFA <- function(data, K){
  data_t <- lapply(data, t)
  MOFAobject <- createMOFAobject(data_t)
  TrainOptions <- getDefaultTrainOptions()
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  ModelOptions$numFactors <- 20
  DataOptions <- getDefaultDataOptions()
  TrainOptions$DropFactorThreshold <- 0.1
  n_inits <- 4
  MOFAlist <- lapply(1:n_inits, function(it) {

    TrainOptions$seed <- 2018 + it

    MOFAobject <- prepareMOFA(
      MOFAobject,
      DataOptions = DataOptions,
      ModelOptions = ModelOptions,
      TrainOptions = TrainOptions
    )
    out <- tryCatch(
      {
        runMOFA(MOFAobject)
      },
      error=function(cond) {
        message(paste("Shut down all components, no structure found in the data."))
        # Choose a return value in case of error
        return(NA)
      }
    )
    return(out)
  })
  ## all is na
  if(sum(sapply(MOFAlist, is.na))==4){
    res <- NA
  }else{
    if(sum(sapply(MOFAlist, is.na))>=1){
      MOFAlist <- MOFAlist[-which(is.na(MOFAlist))]
    }
    MOFAobject <- selectModel(MOFAlist, plotit = FALSE)
    clust.MOFA <- MOFAobject@Expectations$Z%>%
      dist %>% hclust(method="ward.D2") %>% cutree(k=K)
    res <- list(clust=clust.MOFA, fit = MOFAobject)
  }
  return(res)

}
>>>>>>> origin/review
