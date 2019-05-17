#' simulateY
#' @title Simulate various type of data
#'
#' @param nclust A numeric value specifying the number of cluster (integer)
#' @param n_byClust A list or numeric value specifying the number of patients by cluster (if only one value is specify groups contain the same number of patients)
#' @param J A numeric value specifying the number of total biomarkers
#' @param prop A numeric value specifying the proportion of informative biomarkers for each cluster
#' @param noise A numeric value specifying the noise for back ground distribution
#' @param flavor  A character value specifying the distribution of simulations
#'   Defaults to \code{"normal"}. See Details
#' @param params A list of the parameters for flavor (See Details)
#'
#' @details The observations are simulated as follows according to the
#'   value of argument 'flavor':
#'
#'   If \code{flavor=="normal"} (the default), the observations follow a normal distribution (specify in \code{params} the mean and the standard deviance for each cluster.
#'   \code{params} is a list of vector \code{list(c(mean, sd))}. If the length is one we
#'   consider that parameters are the same for all groups, else the length of the list needs
#'   to be of length \code{nclust}
#'
#'   If \code{flavor=="beta"}, the observations follow a beta distribution (methylation data).
#'
#'   If \code{flavor=="binary"}, the observations follow a bernoulli (0 or 1) distribution (specify in \code{params} the proportion of 1). For example to simulate mutations
#'   \code{params} is a list of vector \code{list(c(p))}. If the length is one we
#'   consider that parameters are the same for all groups, else the length of the list needs
#'   to be of length \code{nclust}
#'
#' @return A list with two components: \describe{ \item{data}{A \code{nclust*n_byclust} x
#'   \code{J} matrix}, \item{positive}{the features that drives clusters} }
#' @examples
#' c <- simulateY(J=100, prop=0.1)
#' heatmap(c$data, scale="none")
#' c_bin <- simulateY(J=100, flavor="binary",
#' params=list(c(p=0.7)), prop=0.1, noise=0.1)
#' heatmap(c_bin$data, scale="none")
#' c_beta <- simulateY(J=400, flavor="beta",
#' params=list(c(mean1=-3,mean2=3, sd1=0.1, sd2=0.2) ), prop=0.3, noise=1)
#' heatmap(c_beta$data[,unlist(c_beta$positive)], scale="none")
#' @importFrom dplyr %>%
#' @importFrom Matrix forceSymmetric
#' @import InterSIM
#' @export
simulateY <- function(nclust=4, n_byClust=20, J, prop=0.01, noise=0.1, flavor=c("normal", "beta", "binary"), params=list(c(mean=1,sd=1))){
  ## Sanity checks
  if (!is.numeric(nclust)) {
    stop("nclust must be an integer")
  }
  if (length(n_byClust)==1) {
    n_byClust <- rep(n_byClust, nclust)## only one value same size for all groups
  }
  if(length(n_byClust)!=nclust){
    stop(sprintf("length of n_byClust must be equal to nclust=%s", nclust))
  }
  n = sum(n_byClust)  # total number of samples

  if (length(prop)==1) {
    prop <- rep(prop, nclust)## only one value same prop for all groups
  }
  if(length(prop)!=nclust){
    stop(sprintf("length of prop must be equal to nclust=%s", nclust))
  }

  if (length(params)==1) {
    params <- replicate(nclust,params, simplify = TRUE)## only one value same parameters for all groups
  }
  if(length(params)!=nclust){
    stop(sprintf("number of row of params must be equal to nclust=%s", nclust))
  }
  flavor <- match.arg(flavor)

  true.clusters = mapply(function (clust, n){
    c(rep(sprintf("C%s", clust),n))
  }, 1:nclust, n_byClust, SIMPLIFY = FALSE) %>% unlist

  flavor <- match.arg(flavor)

  if(flavor=="beta"){
    ep <- c( "mean1", "mean2", "sd1",   "sd2"  )
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse="','"))
      stop("Argument 'params' should for normal contain columns named ", str)
    }
  }
  if(flavor=="normal"){
    ep <- c( "mean", "sd")
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse="','"))
      stop("Argument 'params' for beta should contain columns named ", str)
    }
  }
  if(flavor=="binary"){
    ep <- c("p")
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse="','"))
      stop("Argument 'params' for binary should contain columns named ", str)
    }
  }
  data_sim = matrix(rnorm(n*J, mean=0, sd=noise),nrow=n, ncol=J) ## Nul matrix +noise
  simulate_data <- switch(flavor,
                          "normal"=simulate_norm,
                          "beta"=simulate_beta,
                          "binary"=simulate_binary)


  data_sim <- mapply(simulate_data,
                     n_byClust,
                     round(prop*J),
                     params,
                     J, SIMPLIFY = FALSE
  )
  data <- do.call(rbind, lapply (data_sim, function (dd) dd$C))

  if(flavor=="binary"){
    data <- data + matrix(rbinom(n*J,1,noise),nrow=n, ncol=J)
    data <- pmin(data,1) ## Max is one
  }
  else if(flavor == "beta"){
    rev.logit <- function(x)  1/(1+exp(-x+ matrix(rnorm(n*J, mean=0, sd=noise),nrow=n, ncol=J)))
    data <- rev.logit(data)
  }else{
    data <- data + matrix(rnorm(n*J, mean=0, sd=noise),nrow=n, ncol=J)


  }
  ## Named columns
  colnames(data) <- sprintf("gene%s", 1:ncol(data))
  positive <- lapply (data_sim, function (dd) colnames(data)[dd$positive])
  return(list(data=data, positive=positive, true.clusters=true.clusters))
}



################# normal distribution ####################

simulate_norm <- function(n,## patient
                          j, ## biomark
                          params, ## parameter distribution
                          J ## all biomarkers
){
  m <- params["mean"]
  sd <- params["sd"]
  c= matrix(rnorm(n*j,mean=m,sd=sd),ncol=j,nrow=n)
  c= cbind(c, matrix(0,nrow=n, ncol=J-j))
  idx <- sample(1:ncol(c))
  positive <- sapply(1:j, function (ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}
################# distribution for methylation ####################

simulate_beta <- function(n,## patient
                          j, ## biomark
                          params, ## parameters
                          J ## all biomarkers
){
  m1 <- params["mean1"]
  sd1 <- params["sd1"]
  m2 <- params["mean2"]
  sd2 <- params["sd2"]
  c_1 <- matrix(rnorm(n*3*j/4,mean=m1,sd=sd1),ncol=3*j/4,nrow=n) ## Hypo (neg mean)
  c_2 <- matrix(rnorm(n*j/4,mean=m2,sd=sd2),ncol=j/4,nrow=n) ## Hyper (pos mean)

  c = cbind(c_1, c_2)
  c = cbind(c, matrix(logit(runif(n*(J-j))), ncol=J-j, nrow=n))
  idx <- sample(1:ncol(c)) ## sample columns
  positive <- sapply(1:j, function(ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}
################# binary distribution ####################

simulate_binary <- function(n,## patient
                            j,## biomark
                            params,
                            J ## all biomarkers
){
  p <- params["p"]
  c=matrix(rbinom(n*j,1,p),ncol=j,nrow=n)
  c= cbind(c, matrix(0,nrow=n, ncol=J-j))
  idx <- sample(1:ncol(c))
  positive <- sapply(1:j, function (ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}

logit <- function (x) {
  log(x) - log(1 - x)
}

