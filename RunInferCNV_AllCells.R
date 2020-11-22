# Load function
source('Utils.R')

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

library(infercnv)

data <- readRDS(paste0(outDir, 'data_all.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

data <- data[rownames(data)%in%geneLocs[,1],]
sampleAnnotation <- data@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(data@assays$RNA@counts),
                                     annotations_file = sampleAnnotation,
                                     gene_order_file = geneLocs,
                                     ref_group_names = ref_groups)

# Make an output directory for this run
inferCNVOut <- paste0(outDir, 'InferCNV_allCells_', date)
dir.create(inferCNVOut)

# Auto detect the number of cpus
# Otherwise it defaults to 4 and we only utilize 25% of our cloud server XD
nCores <- parallel::detectCores()
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = inferCNVOut,
                              cluster_by_groups = T,
                              HMM = T,
                              num_threads = nCores)
