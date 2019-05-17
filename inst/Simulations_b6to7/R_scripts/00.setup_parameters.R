### This script setups the parameters for the simulation
library(ggplot2)
library("tibble")
library("dplyr")
library(CrIMMix)
library("tidyverse")

## Parallel with future
pathDat <- R.utils::Arguments$getWritablePath("Data_sim_20181012")
pathMeth <- R.utils::Arguments$getWritablePath("Data_Results_20181012")

means <- list(rep(2,2), rep(2, 3))
sds <- list(rep(1,2), rep(1,3))
params_beta <- list(c(mean1=-2, mean2=2, sd1=0.5, sd2=0.5))
S <- 50
nclust=c(2, 3)
n_byClust=list(c(25,25), c(25,25,25))
force=FALSE

grid.param <- list(noiseD1=c(0.2),
                   noiseD2=c(0.1)/10,
                   noiseD3=c(0.1)*3)
props <- c(0.005, 0.01, 0.02)

nnD1 <- grid.param$noiseD1
nnD2 <- grid.param$noiseD2
nnD3 <- grid.param$noiseD3
nbCPU=15L
