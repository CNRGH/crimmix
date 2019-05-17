cv <- readRDS("inst/extdata/TCGA/Liver/cv_liver.rds") #( Run on bioinfo)
library(CrIMMix)
library(dplyr)
library(ggplot2)
idx_dat <- 1:4
count_var.df <- do.call(rbind, lapply(names(cv), function (mm){
  cv_s<- cv[[mm]]

  extract_pos_s<- sapply(cv_s, function (cv_ss){
    fit <- cv_ss$fit
    extractPos(fit, mm)
  })
  count_var <- apply(extract_pos_s, 1, function (ll){
    ll%>% unlist %>% table %>% sort(decreasing=TRUE)/length(cv_s) %>% as.numeric
  })
  count_var <- lapply(count_var, as.numeric)
  n_dat <- sapply(count_var, length)
  count_var_dat <- do.call(c, count_var) %>% as.data.frame
  colnames(count_var_dat) <- "Frequency"
  count_var_dat$dataset <- sprintf("dataset %s", rep(idx_dat  , times=n_dat))
  count_var_dat$method <- mm
  count_var_dat
}))

library(ggplot2)
color <- RColorBrewer::brewer.pal(9, "Spectral")[6:9]
count_var.df$method <- factor(count_var.df$method, levels=c("SGCCA", "Mocluster"))
p <- ggplot(count_var.df, aes(x=method, y=Frequency, fill=method)) +
  geom_violin()+ scale_fill_manual(values=color)+facet_wrap(.~dataset, ncol=2)+theme_bw()
p
