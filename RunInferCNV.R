# Load function
source('Utils.R')

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

library(infercnv)

data <- readRDS(paste0(outDir, 'dataFiltered.RDS'))
geneLocs <- read.table(paste0(dataDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

data.compact <- data.compact[rownames(data.compact)%in%geneLocs[,1],]
sampleAnnotation <- data.compact@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               sampleAnnotation$cell.type)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(data.compact@assays$RNA@counts),
                                     annotations_file = sampleAnnotation,
                                     gene_order_file = geneLocs,
                                     ref_group_names = c('Healthy'))

# Make an output directory for this run
inferCNVOut <- paste0(outDir, 'InferCNV_', date)
dir.create(inferCNVOut)

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = inferCNVOut,
                              cluster_by_groups = T,
                              HMM = T)
