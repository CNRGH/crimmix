#' roc_eval computation
#'
#' @param truth a list of biomarkers for each data set
#' @param fit the result of the method
#' @param method method used
#'
#' @return roc_eval (TPR and FPR) for each data set
#' @export
#' @importFrom stringr str_extract
roc_eval <- function(truth, fit, method){

  doROC <- switch(method,
                  "iCluster" = doROC_iCluster,
                  "MOFA" = doROC_MOFA,
                  "SGCCA" =  doROC_SGCCA,
                  "RGCCA" =  doROC_RGCCA,
                  "NMF" =  doROC_intNMF,
                  "MCIA" =  doROC_MCIA,
                  "Mocluster" = doROC_Moa, 
                  "CIMLR" = doROC_CIMLR
  )
  res <- doROC(truth,fit)
  return(res)
}

doROC_CIMLR <- function(truth, fit){
  regexp <- "[[:digit:]]+"
  selectVars_1 <- fit$selectfeatures$names[order(fit$selectfeatures$pval)]
  k_grid <- stringr::str_extract(pattern="_dat*.",selectVars_1) %>% unlist %>% unique %>% sort()
  test_o <- lapply(k_grid, function (kk){
    idx <- grep(kk, selectVars_1)
    gsub(kk, "", selectVars_1[idx])
  }) %>% lapply(function (ll) ll %>% str_extract(pattern=regexp))
  
  J <- sapply(test_o, length)
  
  denom_tp <- sapply(truth, length)
  TPR_list <-lapply(1:length(test_o), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })
  
  FPR_list <- lapply(1:length(test_o), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  
  return(list(TPR = TPR_list, FPR = FPR_list))
  
}

doROC_iCluster <- function(truth, fit){
  a <- lapply(1:length(fit$beta), function(ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    names(rowsum) <- 1:length(rowsum)
    upper=quantile(rowsum,prob=0.85)
    sigfeatures=names(which(rowsum>upper))
  })
  regexp <- "[[:digit:]]+"
  truth <- lapply(truth, str_extract, pattern=regexp)
  # process string
  test_o <- lapply(1:length(fit$beta), function (ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    names(rowsum) <- 1:length(rowsum)
    names(sort(rowsum[a[[ii]]], decreasing = TRUE))
  })
  J <- sapply(fit$beta, nrow)
  
  denom_tp <- sapply(truth, length)
  TPR_list <-lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]] %>% as.numeric) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })

  return(list(TPR = TPR_list, FPR = FPR_list))
}

doROC_Moa <- function(truth, fit){
  K <- fit@data %>% length
  a <- fit@loading
  J <- sapply(fit@data, nrow)
  selectVars_1 <- which(a %>% rowSums !=0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  regexp <- "[[:digit:]]+"
  
  test <- rowSums(abs(a))
  idx <- which(test!=0)
  test_o <- sort(test[idx], decreasing = TRUE) %>% names
  test_o <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), test_o)
    gsub(sprintf("_dat%s", kk), "", test_o[idx]) %>% str_extract(pattern=regexp)
  })
  denom_tp <- sapply(truth, length)
  
  TPR_list <- lapply(1:K, function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:K, function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}

doROC_SGCCA <- function(truth, fit){
  a <- fit$a
  J <- sapply(a, nrow)
  selectVars <- lapply(a, function(aa) which(rowSums(aa) != 0) %>% names)
  test <- lapply(a, function(aa) rowSums(abs(aa)))
  regexp <- "[[:digit:]]+"
  denom_tp <- sapply(truth, length)
  
  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names %>% str_extract(pattern=regexp)
  })

  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}

doROC_MOFA <- function(truth, fit){
  a <- fit@Expectations$W

  J <- sapply(a, nrow)
  regexp <- "[[:digit:]]+"
  
  selectVars <- lapply(a, function(aa)
    aa %>% apply(2,FUN = function(x) which(abs(x)>1e-2)) %>% unlist %>% rownames
    )
  test <- lapply(a, function(aa) colSums( abs(aa)))
  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names %>% str_extract(pattern=regexp)
  })
  denom_tp <- sapply(truth, length)
  
  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}


doROC_RGCCA <- function(truth, fit){
  a <- fit$a
  J <- sapply(a, nrow)
  regexp <- "[[:digit:]]+"
  
  selectVars <- lapply(a, function(aa) which(rowSums(aa) != 0) %>% names)
  test <- lapply(a, function(aa) rowSums(abs(aa)))
  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names %>% str_extract(pattern=regexp)
  })
  denom_tp <- sapply(truth, length)
  
  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}


doROC_MCIA<- function(truth, fit){
  coas <- fit$coa
  a <- lapply(coas, function (cc) cc$li)
  J <- sapply(a, nrow)
  regexp <- "[[:digit:]]+"
  
  test <- lapply(a, function(aa) rowSums(abs(aa)))
  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names %>% str_extract(pattern=regexp)
  })
  denom_tp <- sapply(truth, length)
  
  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}

doROC_intNMF <- function(truth, fit){
  a <- fit$H
  J <- sapply(a, ncol)
  regexp <- "[[:digit:]]+"
  
  test <- lapply(a, function(aa) colSums((aa)))
  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names %>% str_extract(pattern=regexp)
  })
  denom_tp <- sapply(truth, length)
  
  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}


#' Plot ROC curve
#'
#' @param roc_eval output of \code{roc_eval} function
#' @export
#' @importFrom ggplot2 ggplot  geom_line scale_x_continuous scale_y_continuous theme_bw
plot_roc_eval <- function(roc_eval){
  ROC_df <- do.call(rbind, lapply(1:length(roc_eval$FPR), function (ii){
    fpr <- roc_eval$FPR[[ii]]
    tpr <- roc_eval$TPR[[ii]]
    return(data.frame(tpr=tpr,fpr=fpr, data=sprintf("data%s", ii)))
  }))
  ggplot(ROC_df, aes(x=fpr, y=tpr, color=data, type=data))+geom_line()+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+theme_bw()

}



