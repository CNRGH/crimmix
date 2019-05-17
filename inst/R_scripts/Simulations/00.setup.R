### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### This script setups the parameters for the simulation benchmarks of article 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
library(CrIMMix)
library(dplyr)
## Define path to save results
pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
pathMeth <- R.utils::Arguments$getWritablePath("Data_Results_20181012")
## Number of simulations
S <- 50

## Number of individuals by clust
nclust <- c(4, 4, 4, 4, 4, 2, 3, 4)
n_by_Clust <- list(B1 = c(10,20,5,25), B2 = c(10,20,5,25), B3 = c(10,20,5,25), B4 = c(10,20,5,25), B5 = c(10,20,5,25), B6 = c(25,25), B7 = c(25,25,25), B8=c(25,25,25,25))

## Length of datasets
J_gaussian <- rep(1000, length(nclust))
J_beta <- rep(5000, length(nclust))
J_binary <- rep(500, length(nclust))

## Define means and sds for gaussian dataset
means_grid <- list(B1 = rep(2,4), B2 = rep(2,4), B3 = rep(2,4), B4 = rep(2,4), B5 = rep(2,4), B6 = rep(2,2), B7 = rep(2,3), B8 = rep(2,4))
sds_grid <- list(B1 = rep(1,4), B2 = rep(1,4), B3 = rep(1,4), B4 = rep(1,4), B5 = rep(1,4), B6 = rep(1,2), B7 = rep(1,3), B8 = rep(1,4))
params_gaussian <- lapply(1:length(means_grid),function(ii) mapply(function(m, sd) return(c(mean = m, sd = sd)), means_grid[[ii]], sds_grid[[ii]], SIMPLIFY = FALSE))

## define parameters for beta 
params_beta <- lapply(1:8, function(ii) c(mean1 = -2, mean2 = 2, sd1 = 0.5, sd2 = 0.5))
names(params_beta) <- sprintf("B%s", 1:8)

## define parameters for binary
params_bin <- lapply(1:8, function(ii)  list(c(p = 0.6)))
names(params_beta) <- sprintf("B%s", 1:8)


## Define noise
B1_props = c(0.005, 0.01, 0.02)
multiple <- c(1, 1, 1, 1, 2, 1, 1,1)
props <- lapply(multiple, function(mm) B1_props*mm)
names(props) <- sprintf("B%s", 1:8)

noiseD1 <- c(0.2, 0.5, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2)
noiseD2 <- c(0.1, 0.1, 0.5, 0.2, 0.1, 0.1, 0.1, 0.1)/10
noiseD3 <- c(0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1)*3



