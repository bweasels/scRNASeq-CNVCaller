library(infercnv)

# Load function
source('Utils.R')

# Try to run CNVcaller on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

### TODO: add functions in correct order after creating inferCNV object