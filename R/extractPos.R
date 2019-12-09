<<<<<<< HEAD
#' Extract positive biomarkers
#'
#' @param fit the result of the method
#' @param method method used
#'
#' @return selected Variables
#' @export
#'
extractPos <- function(fit, method){

  extract <- switch(method,
                 "iCluster"=extract_iCluster,
                 "MOFA"=extract_MOFA,
                 "SGCCA"= extract_SGCCA,
                 "Mocluster"=extract_Moa,
  )
  res <- extract(fit)
  return(res)
}


extract_iCluster <- function(fit){
  selectVars <- lapply(1:length(fit$beta), function(ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    upper=quantile(rowsum,prob=0.75)
    sigfeatures=which(rowsum>upper)
  })
  return(selectVars)
}

extract_Moa <- function(fit){
  K <- fit@data %>% length
  a <- fit@loading
  selectVars_1 <- which(a %>% rowSums !=0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  return(selectVars)
}

extract_SGCCA <- function(fit){
  a <- fit$a
  selectVars <- lapply(a, function(aa) which(rowSums(aa) != 0) %>% names)
  return(selectVars)
}

extract_MOFA <- function( fit){
  a <- fit@Expectations$W
  selectVars <- lapply(a, function(aa)
    aa %>% apply(2,FUN = function(x) which(abs(x)>1e-2)) %>% unlist %>% unique %>% names
  )
  return(selectVars)
}
||||||| merged common ancestors
=======
#' Extract positive biomarkers
#'
#' @param fit the result of the method
#' @param method method used
#'
#' @return selected Variables
#' @export
#'
extractPos <- function(fit, method){

  extract <- switch(method,
                 "iCluster"=extract_iCluster,
                 "MOFA"=extract_MOFA,
                 "SGCCA"= extract_SGCCA,
                 "Mocluster"=extract_Moa,
                 "CIMLR"=extract_CIMLR
  )
  res <- extract(fit)
  return(res)
}

extract_CIMLR <- function(fit){
    selectVars_1 <- fit$selectfeatures$names[fit$selectfeatures$pval<1e-5]
    k_grid <- stringr::str_extract(pattern="_dat*.",selectVars_1) %>% unlist %>% unique %>% sort()
    selectVars <- lapply(k_grid, function (kk){
      idx <- grep(kk, selectVars_1)
      gsub(kk, "", selectVars_1[idx])
    })
  return(selectVars)
}

extract_iCluster <- function(fit){
  selectVars <- lapply(1:length(fit$beta), function(ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    upper=quantile(rowsum,prob=0.75)
    sigfeatures=which(rowsum>upper)
  })
  return(selectVars)
}

extract_Moa <- function(fit){
  K <- fit@data %>% length
  a <- fit@loading
  selectVars_1 <- which(a %>% rowSums !=0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  return(selectVars)
}

extract_SGCCA <- function(fit){
  a <- fit$a
  selectVars <- lapply(a, function(aa) which(rowSums(aa) != 0) %>% names)
  return(selectVars)
}

extract_MOFA <- function( fit){
  a <- fit@Expectations$W
  selectVars <- lapply(a, function(aa)
    aa %>% apply(2,FUN = function(x) which(abs(x)>1e-2)) %>% unlist %>% unique %>% names
  )
  return(selectVars)
}
>>>>>>> origin/review
