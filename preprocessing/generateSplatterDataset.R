## Generate our test dataset with Splatter

# Use source for homebrew scripts
source('Utils.R')

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

data <- readRDS(paste0(outDir, 'data_all.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)
# Trimming dataset to genes available in GeneLocs (from 19184 to 15252 genes)
data <- data[rownames(data)[(rownames(data) %in% rownames(geneLocs))]]

# Get essential stats for splatter simulation
if (!file.exists(paste0(outDir, 'InitialParams_2020-11-24.RDS'))) {
  params <- splatEstimate(as.matrix(data@assays$RNA@counts))
  saveRDS(params, paste0('InitialParams_', date, '.RDS'))
} else {
  params <- readRDS('InitialParams_2020-11-24.RDS')
}

# params: ngenes = 15252, ncells = 5278
# Need to figure out parameters to mock cnv variations
sim <- splatSimulate(params)
saveRDS(sim, paste0(outDir, 'data_sim.RDS'))
