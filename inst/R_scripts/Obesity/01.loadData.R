library(dplyr)
library(CrIMMixdata)
data(obesity)

idx_AER <- grep("AER", rownames(obesity[[1]]))
idx_PRT <- grep("PRT", rownames(obesity[[1]]))
idx_POST <- grep("POST", rownames(obesity[[1]]))
## Filtering
### Methylome
meth <- obesity[[1]]
hypo <- which(colSums(meth<=0.1)== nrow(meth))
hyper <- which(colSums(meth>=0.9)== nrow(meth))
meth <- meth[, -c(hypo, hyper)]


obesity_filter <- obesity
obesity_filter[[1]] <- meth
idx_pat <- seq(from=1, by =2, length=13)
obesity_filter_scale <- lapply(obesity_filter, scale)

new_dat <- lapply(obesity_filter_scale, function (jj){
    n_d <- sapply(idx_pat, function(ii){
      jj[ii,]-jj[ii+1,]
}) %>% t
    rownames(n_d) <- rownames(jj[idx_pat,])
    return(n_d)
  })
str(new_dat)

