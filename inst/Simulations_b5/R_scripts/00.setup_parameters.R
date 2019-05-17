### This script setups the parameters for the simulation
library(ggplot2)
library("tibble")
library("dplyr")
library(CrIMMix)
library("tidyverse")

## Parallel with future
pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
pathMeth <- R.utils::Arguments$getWritablePath("Data_Results_20181012")

means <- c(2,2,2,2)
sds <- c(1,1,1,1)
params <- mapply(function (m, sd) return(c(mean=m, sd=sd)), means, sds, SIMPLIFY=FALSE)
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 50
nclust=4
n_byClust=c(10,20,5,25)
force=FALSE

grid.param <- list(noiseD1=c(0.2),
                   noiseD2=c(0.1)/10,
                   noiseD3=c(0.1)*3)
props <- c(0.005, 0.01, 0.02)*2

nnD1 <- grid.param$noiseD1
nnD2 <- grid.param$noiseD2
nnD3 <- grid.param$noiseD3
nbCPU=15L
