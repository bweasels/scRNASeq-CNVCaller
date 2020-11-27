# Load function
source('Utils.R')

# Prep for HoneyBADGER on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

date <- Sys.Date()

library(HoneyBADGER)
library(Seurat)
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)
library(parallel)

# list data paths
data.paths <- list.files(path = outDir,
                         pattern = '_UMM.*.RDS',
                         full.names = T)
sampleNames <- gsub('.*Melanoma_(.*).RDS$', '\\1', data.paths)
inputs <- list(list(data.paths[1], sampleNames[1], plottingDir),
               list(data.paths[2], sampleNames[2], plottingDir))

# Running in parallel on both samps b/c USE ALL THE CORES!!
RunHoneyBADGER <- function(params){
  data.path <- params[[1]]
  sample.name <- params[[2]]
  plot.dir <- params[[3]]
  
  # Load data, and the start stop locations of each transcript
  data <- readRDS(data.path)
  
  # Trim tx in data to match tx I was able to find (~500 genes different - all ribosomal varants)
  sampleAnnotation <- data@meta.data
  sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                                 cellType = sampleAnnotation$cell.type)
  
  # Get the unique cell types and remove HCC so we have a list of reference cell types
  ref_groups <- unique(sampleAnnotation$cellType)
  ref_groups <- ref_groups[!grepl('Melanoma', ref_groups)]
  
  Mel <- data[,data$cell.type=='Melanoma']
  Mel <- as.matrix(Mel@assays$RNA@counts)
  
  ref <- data[,data$cell.type!='Melanoma']
  ref <- as.matrix(ref@assays$RNA@counts)
  
  # Get the gene positions
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = 'hsapiens_gene_ensembl')
  
  # make honeybadger object
  hb <- new('HoneyBADGER', name=sample.name)
  hb$setGexpMats(gexp.sc.init = Mel, 
                 gexp.ref.init = ref, 
                 mart.obj = mart.obj, 
                 filter = F, 
                 scale = T, 
                 verbose = T)
  
  # Model variance and necessary exp deviance to ID CNV
  hb$setMvFit(verbose = T)
  hb$setGexpDev(verbose=T)
  hb$calcGexpCnvBoundaries(init=T,
                           verbose = T)
  
  # Output the identified CNVs
  print(range(hb$genes[unlist(hb$bound.genes.final)]))
  
  # Retest to ID posterior probability & get final results
  hb$retestIdentifiedCnvs(retestBoundGenes = T,
                          retestBoundSnps = F,
                          verbose = T)
  
  results <- hb$summarizeResults(geneBased = T,
                                 alleleBased = F)
  
  # Make result plot
  pdf(paste0(plottingDir, 'HoneyBADGERRes_', sample.name, '_', date, '.pdf'))
  trees <- hb$visualizeResults(geneBased = T,
                               alleleBased = F,
                               details = T,
                               margins = c(25, 15))
  
  hc <- trees$hc
  order <- hc$labels[hc$order]
  
  hb$plotGexpProfile()
  hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['amp']])
  hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['del']])
  dev.off()
}

mclapply(inputs, RunHoneyBADGER)
