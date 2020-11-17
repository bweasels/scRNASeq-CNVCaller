## Generate our test datasets with Splatter

# Use source for homebrew scripts
source('SetDirectory.R')

# @Charlotte Please enter your own ifelse statement in this function
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

#Load libraries and set date
library(splatter)
library(Seurat)

date <- Sys.Date()

#set.seed locks our random number indexer so we get repeatable results
set.seed(10000)

data <- readRDS(paste0(outDir, 'dataFiltered_2020-11-15.RDS'))

# Get essential stats for splatter simulation
# 11K Cells took ~1.5 hours on my laptop
if(!file.exists(outDir, 'InitialParams_2020-11-16.RDS')){
  params <- splatEstimate(as.matrix(data@assays$RNA@counts))
  saveRDS(params, paste0('InitialParams_', date, '.RDS'))
}else{
  params <- readRDS('InitialParams_2020-11-16.RDS')
}

# Need to figure out parameters to mock cnv variations
sim <- splatSimulate(params)