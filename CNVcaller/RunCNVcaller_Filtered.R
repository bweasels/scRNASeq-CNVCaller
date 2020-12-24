library(infercnv)

source('CNVcaller/CNVcaller_Algo')
source('Utils.R')

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

# Load data, and the start stop locations of each transcript
data <- readRDS(paste0(outDir, 'dataFiltered.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

# Trim tx in data to match tx I was able to find (~500 genes different - all ribosomal varants)
data <- data[rownames(data)[(rownames(data) %in% rownames(geneLocs))]]
sampleAnnotation <- data@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

# Make the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = as.matrix(data@assays$RNA@counts),
                                     annotations_file = sampleAnnotation,
                                     gene_order_file = geneLocs,
                                     ref_group_names = ref_groups)

# Make an output directory for this run
CNVcallerOut <- paste0(outDir, 'CNVcaller_Filtered_test_', date)
dir.create(CNVcallerOut)

# Auto detect the number of cpus
# Otherwise it defaults to 4 and we only utilize 25% of our cloud server XD
nCores <- parallel::detectCores()

# Run infercnv
#infercnv_obj <- refactored_run(infercnv_obj,
infercnv_obj <- run_no_smoothing(infercnv_obj,
                                 cutoff = 0.1,
                                 out_dir = CNVcallerOut,
                                 cluster_by_groups = T,
                                 HMM = T,
                                 num_threads = nCores,
                                 denoise = T)